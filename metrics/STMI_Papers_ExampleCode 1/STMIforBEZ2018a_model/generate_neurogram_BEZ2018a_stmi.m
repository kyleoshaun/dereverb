function [neurogram_stmi,t_stmi,CFs] = generate_neurogram_BEZ2018a_stmi(stim,Fs_stim,species,ag_fs,ag_dbloss)
%GENERATE_NEUROGRAM_BEZ2018_STMI Generate a neurogram with CFs and time
% steps suitable for the STMI speech-intelligibility metric, using the
% auditory periphery model of Bruce, Erfani & Zilany (Hear. Res. 2018),
% including an option to use a humanized version of the model.
%
% [NEUROGRAM_STMI,T_STMI,CFS] = GENERATE_NEUROGRAM_BEZ2018_STMI(STIM,FS_STIM,SPECIES,AG_FS,AG_DBLOSS)
%
% NEUROGRAM_STMI is the stmi neurogram matrix
% T_STMI is the array of times (in second) for the stmi neurogram
% CFS is the array of CFs in Hz for both of the neurograms
% STIM is the time-domain stimulus waveform
% FS_STIM is the sampling frequency of the stimulus
% SPECIES is the model species: 1 = cat, 2 = human
% AG_FS and DBLOSS are the audiogram frequencies (in Hz) and corresponding
%  hearing thresholds (in dB HL)
%
% Ian C. Bruce <ibruce@ieee.org> © 2024


% model fiber parameters
CFs   = 440 * 2.^ (((0:127)-31)/24);  % Logarithmically-spaced CFs in Hz
numcfs = length(CFs); % Number of CFs to use in the neurogram

% cohcs  = ones(1,numcfs);  % normal ohc function
% cihcs  = ones(1,numcfs);  % normal ihc function

dbloss = interp1(ag_fs,ag_dbloss,CFs,'linear','extrap');

% mixed loss
[cohcs,cihcs,OHC_Loss]=fitaudiogram2(CFs,dbloss,species);

% OHC loss
% [cohcs,cihcs,OHC_Loss]=fitaudiogram(CFs,dbloss,dbloss);

% IHC loss
% [cohcs,cihcs,OHC_Loss]=fitaudiogram(CFs,dbloss,zeros(size(CFs)));

numsponts_healthy = [10 10 30]; % Number of low-spont, medium-spont, and high-spont fibers at each CF in a healthy AN

if exist('ANpopulation.mat','file')
    load('ANpopulation.mat');
    disp('Loading existing population of AN fibers saved in ANpopulation.mat')
    if (size(sponts.LS,2)<numsponts_healthy(1))||(size(sponts.MS,2)<numsponts_healthy(2))||(size(sponts.HS,2)<numsponts_healthy(3))||(size(sponts.HS,1)<numcfs||~exist('tabss','var'))
        disp('Saved population of AN fibers in ANpopulation.mat is too small - generating a new population');
        [sponts,tabss,trels] = generateANpopulation(numcfs,numsponts_healthy);
    end
else
    [sponts,tabss,trels] = generateANpopulation(numcfs,numsponts_healthy);
    disp('Generating population of AN fibers, saved in ANpopulation.mat')
end

implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
noiseType = 1;  % 0 for fixed fGn (1 for variable fGn)
expliketype = 1; % 1 for shifted softplus (preferred); 0 for no expontential-like function; 2 for shifted exponential; 3 for shifted Boltmann

% stimulus parameters
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
stim100k = resample(stim,Fs,Fs_stim).';
T  = length(stim100k)/Fs;  % stimulus duration in seconds

% PSTH parameters
nrep = 1;
psthbinwidth_stmi = 8e-3; % stmi neurogram time binwidth in seconds
windur_stmi=2;
smw_stmi = rectwin(windur_stmi); % Smoothing window for the stmi neurogram

% Calculate a simulation duration of approximately 1.2 times the stimulus
% duration, to allow for some recovery to the spontaneous rate after the
% stimulus, making sure that the length is an integer multiple of the
% stmi neurogram time binwidth
simdur = ceil(T*1.2/psthbinwidth_stmi)*psthbinwidth_stmi;

pin = stim100k(:).';

clear stim100k

for CFlp = 1:numcfs
    
    
    CF = CFs(CFlp);
    cohc = cohcs(CFlp);
    cihc = cihcs(CFlp);
    
    numsponts = round([1 1 1].*numsponts_healthy); % Healthy AN
    %     numsponts = round([0.5 0.5 0.5].*numsponts_healthy); % 50% fiber loss of all types
    %     numsponts = round([0 1 1].*numsponts_healthy); % Loss of all LS fibers
    %     numsponts = round([cihc 1 cihc].*numsponts_healthy); % loss of LS and HS fibers proportional to IHC impairment
    
    sponts_concat = [sponts.LS(CFlp,1:numsponts(1)) sponts.MS(CFlp,1:numsponts(2)) sponts.HS(CFlp,1:numsponts(3))];
    tabss_concat = [tabss.LS(CFlp,1:numsponts(1)) tabss.MS(CFlp,1:numsponts(2)) tabss.HS(CFlp,1:numsponts(3))];
    trels_concat = [trels.LS(CFlp,1:numsponts(1)) trels.MS(CFlp,1:numsponts(2)) trels.HS(CFlp,1:numsponts(3))];
    
    vihc = model_IHC_BEZ2018a(pin,CF,nrep,1/Fs,simdur,cohc,cihc,species);
    
    for spontlp = 1:sum(numsponts)
        
        disp(['CFlp = ' int2str(CFlp) '/' int2str(numcfs) '; spontlp = ' int2str(spontlp) '/' int2str(sum(numsponts))])

        % flush the output for the display of the coutput in Octave
        if exist ('OCTAVE_VERSION', 'builtin') ~= 0
            fflush(stdout);
        end
        
        
        spont = sponts_concat(spontlp);
        tabs = tabss_concat(spontlp);
        trel = trels_concat(spontlp);
        
        psth_ft = model_Synapse_BEZ2018a(vihc,CF,nrep,1/Fs,noiseType,implnt,spont,tabs,trel,expliketype);
        psthbins = round(psthbinwidth_stmi*Fs);  % number of psth_ft bins per psth bin
        psth_stmi = sum(reshape(psth_ft,psthbins,length(psth_ft)/psthbins));

        if spontlp == 1
            neurogram_stmi(CFlp,:) = filter(smw_stmi,1,psth_stmi);
        else
            neurogram_stmi(CFlp,:) = neurogram_stmi(CFlp,:)+filter(smw_stmi,1,psth_stmi);
        end
 
    end

end

neurogram_stmi = neurogram_stmi(:,1:windur_stmi/2:end); % 50% overlap in rectangular window
t_stmi = 0:windur_stmi/2*psthbinwidth_stmi:(size(neurogram_stmi,2)-1)*windur_stmi/2*psthbinwidth_stmi; % time vector for the stmi neurogram      
