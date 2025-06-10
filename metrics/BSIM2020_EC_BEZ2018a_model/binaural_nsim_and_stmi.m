function [STMI, NSIM_FT, ec_out] = binaural_nsim_and_stmi(teststim_left, teststim_right, refstim, Fs_stim, stimdb)
% Input data: Pa at the ear drum

enable_plots = true;

%% Add Paths


% EC (BSIM 2020)
addpath('../BSIM_2020')
addpath('../BSIM_2020/SRMRToolbox-master')
addpath('../BSIM_2020/SRMRToolbox-master/libs/')
addpath('../BSIM_2020/SRMRToolbox-master/libs/vad/')
addpath('../BSIM_2020/SRMRToolbox-master/libs/PreProc/')
addpath('../BSIM_2020/SRMRToolbox-master/libs/Gammatonegram/')
addpath('../BSIM_2020/SRMRToolbox-master/libs/Auditory/')

addpath(genpath('../BSIM_2020/Anechoic/'))
addpath('../BSIM_2020/functions/')
addpath('../BSIM_2020/functions/gammatonefilter/')
addpath('../BSIM_2020/functions/ltfat/')


% Model (BEZ2018a)
addpath('../BEZ2018a_model/')

% Other
addpath('../NSIM_Papers_ExampleCode 2/NSIMforBEZ2018a_model/')
addpath('../STMI_Papers_ExampleCode 1/STMIforBEZ2018a_model')


%% Init EC (BSIM 2020)

restore_pwd = pwd;
cd ~/Documents/School/MASc_Thesis/code/metrics/BSIM_2020/
Model_init;
cd(restore_pwd)

%----------------Model parameters-----------------------------------------%
model_param.fs = Fs_stim; % sampling rate 
model_param.fmin = 150; % lowest filter of gammatone filterbank
model_param.fmax = 8500; % highest filter of gammatone filterbank --> Model runs at 100 kHz, but we dont actually need 50kHz BW in the signal right?
model_param.f_target = 500; % specified frequeny that will have a matched filter
model_param.bin_err = 1; % enable(1)/disable(0) binaural processing inaccuracies (important for Monte-Carlo Simulations)
model_param.ERB_factor = 1;% define bandwidth of filters
model_param.long_term = 1; % use long-term model (1) or short-term model (0)
model_param.Filterbank = 'GT';% define filtering: GT is gammatone filter
% Model parameters don't have to be specified. These are also the default values,
% which are used in the model
% Please note that the modified SRMR model, which performs the channel
% selection either considers the whole signal as a single time frame (long
% term or uses the same time constants as the EC processing and BE
% processing if the short time model is used.
%-------------------------------------------------------------------------%

%% Init Model (BEZ2018a)

% Set "reference" audiogram
ag_fs_ref = [0 20e3];
ag_dbloss_ref = [0 0]; % Normal hearing

% Set "test" audiogram
ag_fs_test = [125 250 500 1e3 2e3 4e3 8e3];
%ag_dbloss_test = [0 0 0 0 0 0 0]; % Normal hearing
ag_dbloss_test = [0 0 0 20 40 60 80]; % Example high-freq hearing loss

species = 2; % Human

%% Other Params

% NSIM Params
weights.alpha = 1.0;
weights.beta = 0.0;
weights.gamma = 1.0;
window_type = 0;
win3x3 = ones(3,3);

% STMI Params
rv = 2.^(1:0.5:5);  % cortical (temporal) modulation filter rates
sv = 2.^(-2:0.5:3); % cortical (spectral) modulation filter scales

% Parameters for STMI/NSIM Synthesis
% Regression of Fine-Timing NSIM + STMI performed by Rationalized Arcsine Transformation (RAU)
% SI_SYN = rau_b0 + rau_b1*STMI + rau_b2*NSIM_FT + rau_b3*STMI*NSIM_FT;
% rau_b0 = 0;
% rau_b1 = 1.357;
% rau_b2 = 0.444;
% rau_b3 = 0.00598;


%% Setup Input Data

% Outside function: Set dB SPL
% stimdb = 65; % speech level in dB SPL
% SNR = 0; % in dB
% SNR = inf; % in dB; inf -> no background noise

% refstim = refstim/rms(refstim)*20e-6*10^(stimdb/20);


%% Run EC Front-End

% TEMP
%teststim_left = resample(teststim_left,model_param.fs,Fs_stim);
%refstim = resample(refstim,model_param.fs,Fs_stim);
% END

if length(teststim_left) ~= length(teststim_right)
    error("Left and Right Signal lengths must match")
end

mixed_input = [teststim_left teststim_right];
out_struct = BSIM20('RequiredSignal',mixed_input, 'model_params',model_param);%,...
                  %'OptionalSignal',OptionalSignals,'model_params',model_param);

teststim = out_struct.signals.MixsigSyn'; %out_struct.signals.MixsigSyn;
refstim  = [refstim ; 0];

ec_out = teststim; % Output for reference

%% FT NSIM

% SETUP

if exist('parfor','builtin') % check if the Matlab Parallel Computation
                             % Toolbox is installed and use appropriate
                             % function
    generate_neurogram_function = @generate_neurogram_BEZ2018a_parallelized;
    disp('Using parallelized version of neurogram generation function')
else
    generate_neurogram_function = @generate_neurogram_BEZ2018a;
    disp('Using serial version of neurogram generation function')
end

%refstim = 0;
%teststim = 0;

% NSIM CALCULATION

[neurogram_ft_ref,neurogram_mr_ref,neurogram_Sout_ref,t_ft,t_mr,t_Sout,CFs] = generate_neurogram_function(refstim,Fs_stim,species,ag_fs_ref,ag_dbloss_ref);

[neurogram_ft_test,neurogram_mr_test,neurogram_Sout_test,t_ft,t_mr,t_Sout,CFs] = generate_neurogram_function(teststim,Fs_stim,species,ag_fs_test,ag_dbloss_test);

response_endtime = length(refstim)/Fs_stim+10e-3; % Calculate end time of the neural response to the stimulus
[tmp, ind_mr] = min(abs(t_mr-response_endtime)); % Find the corresponding time index in the mr neurogram
[tmp, ind_ft] = min(abs(t_ft-response_endtime)); % Find the corresponding time index in the ft neurogram

% Scaling method used by Hines and Harte (Speech Comm 2010, 2012)
% scl_mr = 255/max(max(neurogram_mr_ref(:,1:ind_mr)));
% scl_ft = 255/max(max(neurogram_ft_ref(:,1:ind_ft)));

% New scaling method developed by M. R. Wirtzfeld (see Wirtzfeld et al., JARO 2017)
scl_mr = 1/50/t_mr(2);
scl_ft = 1/50/t_ft(2);

[NSIM_MR, ssim_mr] = mssim_v2(scl_mr*neurogram_mr_ref(:,1:ind_mr), scl_mr*neurogram_mr_test(:,1:ind_mr), weights, win3x3, window_type );

[NSIM_FT, ssim_ft] = mssim_v2(scl_ft*neurogram_ft_ref(:,1:ind_ft), scl_ft*neurogram_ft_test(:,1:ind_ft), weights, win3x3, window_type );

% NSIM PLOTS

if enable_plots
    ng1=figure;
    set(ng1,'renderer','painters');
    winlen = 256; % Window length for the spectrogram analyses
    
    sp1 = subplot(3,2,1);
    [s,f,t] = specgram([refstim; eps*ones(round(t_mr(end)*Fs_stim)-length(refstim),1)],winlen,Fs_stim,winlen,0.25*winlen);
    imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
    axis xy; axis tight;
    hcb = colorbar;
    set(get(hcb,'ylabel'),'string','SPL')
    caxis([stimdb-80 stimdb])
    ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
    xlabel('Time');
    ylabel('Frequency (kHz)');
    title('Reference Stimulus Spectrogram')
    xl = xlim;
    
    sp2 = subplot(3,2,2);
    [s,f,t] = specgram([teststim; eps*ones(round(t_mr(end)*Fs_stim)-length(teststim),1)],winlen,Fs_stim,winlen,0.25*winlen);
    imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
    axis xy; axis tight;
    hcb = colorbar;
    set(get(hcb,'ylabel'),'string','SPL')
    caxis([stimdb-80 stimdb])
    ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');
    title('Test Stimulus Spectrogram')
    xlim(xl)
    
    sp3=subplot(3,2,3);
    plot_neurogram(t_mr,CFs,neurogram_mr_ref,sp3);
    caxis([0 80])
    title('Mean-rate Neurogram for Reference Stimulus')
    xlim(xl)
    
    sp4=subplot(3,2,4);
    plot_neurogram(t_mr,CFs,neurogram_mr_test,sp4);
    caxis([0 80])
    title(['Mean-rate Neurogram for Test Stimulus: NSIM\_MR = ' num2str(NSIM_MR)])
    xlim(xl)
    
    sp5=subplot(3,2,5);
    plot_neurogram(t_ft,CFs,neurogram_ft_ref,sp5);
    caxis([0 20])
    title('Fine-timing Neurogram for Reference Stimulus')
    xlim(xl)
    
    sp6=subplot(3,2,6);
    plot_neurogram(t_ft,CFs,neurogram_ft_ref,sp6);
    caxis([0 20])
    title(['Fine-timing Neurogram for Test Stimulus: NSIM\_FT = ' num2str(NSIM_FT)])
    xlim(xl)
end


%% STMI

% STMI CALC SETUP

if exist('parfor','builtin') % check if the Matlab Parallel Computation
                             % Toolbox is installed and use appropriate
                             % function
    generate_neurogram_function = @generate_neurogram_BEZ2018a_stmi_parallelized;
    disp('Using parallelized version of neurogram generation function')
else
    generate_neurogram_function = @generate_neurogram_BEZ2018a_stmi;
    disp('Using serial version of neurogram generation function')
end

rv = 2.^(1:0.5:5);  % cortical (temporal) modulation filter rates
sv = 2.^(-2:0.5:3); % cortical (spectral) modulation filter scales

% refstim = 0;
% teststim = 0;

% STMI CALC

[neurogram_stmi_ref,t_stmi,CFs] = generate_neurogram_function(refstim,Fs_stim,species,ag_fs_ref,ag_dbloss_ref);
base_sig_ref = generate_base_signal(refstim);
[neurogram_base_ref,t_stmi,CFs] = generate_neurogram_function(base_sig_ref,Fs_stim,species,ag_fs_ref,ag_dbloss_ref);

[neurogram_stmi_test,t_stmi,CFs] = generate_neurogram_function(teststim,Fs_stim,species,ag_fs_test,ag_dbloss_test);
base_sig_test = generate_base_signal(teststim);
[neurogram_base_test,t_stmi,CFs] = generate_neurogram_function(base_sig_test,Fs_stim,species,ag_fs_test,ag_dbloss_test);

response_endtime = length(refstim)/Fs_stim+10e-3; % Calculate end time of the neural response to the stimulus
[tmp, ind_stmi] = min(abs(t_stmi-response_endtime)); % Find the corresponding time index in the stmi neurogram

% Generate the cortical responses to the reference signal and reference base signal neurograms and
% compute the reference "template" T matrix
a1_stmi_ref = abs(ngram2cortex(neurogram_stmi_ref(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
a1_base_ref = abs(ngram2cortex(neurogram_base_ref(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
T  = max(a1_stmi_ref - a1_base_ref,0);
TT = sum(sum(T(:).^2));

% Generate the cortical responses to the test signal and test base signal neurograms and
% compute the test "noisy" N matrix
a1_stmi_test = abs(ngram2cortex(neurogram_stmi_test(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
a1_base_test = abs(ngram2cortex(neurogram_base_test(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
N  = max(a1_stmi_test - a1_base_test,0);
NN = sum(sum((max(T(:)-N(:),0)).^2));

STMI = 1-NN/TT;

% STMI PLOTS

if enable_plots
    ng1=figure;
    set(ng1,'renderer','painters');
    winlen = 256; % Window length for the spectrogram analyses
    
    sp1 = subplot(2,2,1);
    [s,f,t] = specgram([refstim; eps*ones(round(t_stmi(end)*Fs_stim)-length(refstim),1)],winlen,Fs_stim,winlen,0.25*winlen);
    imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
    axis xy; axis tight;
    hcb = colorbar;
    set(get(hcb,'ylabel'),'string','SPL')
    caxis([stimdb-80 stimdb])
    ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
    xlabel('Time');
    ylabel('Frequency (kHz)');
    title('Reference Stimulus Spectrogram')
    xl = xlim;

    sp2 = subplot(2,2,2);
    [s,f,t] = specgram([teststim; eps*ones(round(t_stmi(end)*Fs_stim)-length(teststim),1)],winlen,Fs_stim,winlen,0.25*winlen);
    imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
    axis xy; axis tight;
    hcb = colorbar;
    set(get(hcb,'ylabel'),'string','SPL')
    caxis([stimdb-80 stimdb])
    ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');
    title('Test Stimulus Spectrogram')
    xlim(xl)
    
    sp3=subplot(2,2,3);
    plot_neurogram(t_stmi,CFs,neurogram_stmi_ref,sp3);
    caxis([0 200])
    title('STMI Neurogram for Reference Stimulus')
    xlim(xl)
    
    sp4=subplot(2,2,4);
    plot_neurogram(t_stmi,CFs,neurogram_stmi_test,sp4);
    caxis([0 200])
    title(['STMI Neurogram for Test Stimulus: STMI = ' num2str(STMI)])
    xlim(xl)
end


%% Synthesize NSIM/STMI
% Regression of Fine-Timing NSIM + STMI performed by Rationalized Arcsine Transformation (RAU)
% Wirtzfeld et al (2017) 
% https://pmc.ncbi.nlm.nih.gov/articles/PMC5612921/

%SI_SYN = rau_b0 + rau_b1*STMI + rau_b2*NSIM_FT + rau_b3*STMI*NSIM_FT;

end