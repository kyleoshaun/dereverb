clear;

% parameter to explore

expliketype = 1; % sets type of expontential-like function in IHC -> synapse mapping: 1 for shifted softplus (preferred); 0 for no expontential-like function; 2 for shifted exponential; 3 for shifted Boltmann
spontstr = 'HSR'; % sets fiber's spont rate type: 'HSR' for high spont; 'MSR' for medium spont; 'LSR' for low spont
stimdbs = -20:20:60; % stimulus levels in dB SPL
numdbs = length(stimdbs); % number of stimulation levels

%

switch expliketype
    case 1
        explikestr = 'softplus';
    case 0
        explikestr = 'no exp-like func';
    case 2
        explikestr = 'true exp';
    case 3
        explikestr = 'Boltzmann';
    otherwise
        error('Invalid expliketype value')
end

switch spontstr
    case 'HSR'
        % Example high spont rate fiber
        F0 = 600; % stimulus frequency in Hz
        CF = 591;  % fiber characteristic frequency in Hz
        spont = 37; % spont rate in spikes/s

    case 'MSR'

        % Example medium spont rate fiber
        F0 = 700; % stimulus frequency in Hz
        CF = 700;  % fiber characteristic frequency in Hz
        spont = 12; % spont rate in spikes/s

    case 'LSR'

        % Example low spont rate fiber
        F0 = 550; % stimulus frequency in Hz
        CF = 550;  % fiber characteristic frequency in Hz
        spont = 0.1; % spont rate in spikes/s

    otherwise

        error('Unknown spont rate class')

end

% further stimulus parameters

Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
stimdur = 1; % stimulus duration in seconds
stimpts = round(stimdur*Fs); % stimulus duration in samples
restdur = 0.5; % inter-stimulus silence duration in seconds
restpts = round(restdur*Fs); % inter-stimulus silence duration in samples
reppts = stimpts+restpts; % total repeat time in samples
sim_dur = stimdur + restdur; % total simulation duration in seconds
nrep_stim = 400; % number of stimulus repetitions

rt = 5e-3; % rise/fall time in seconds
irpts = round(rt*Fs); % rise/fall time in samples
pin = zeros(1,stimpts); % create stimulus array for one repetition
t = 0:1/Fs:stimdur-1/Fs; % time array for one repetition in seconds
T0 = 1/F0; % stimulus period in seconds

% further AN model parameters
tabs = 0.7e-3; % absolute refractory period in seconds
trel = 0.7e-3; % relative refractory period in seconds
cohc  = 1.0;    % ohc function; 1 = normal
cihc  = 1.0;    % ihc function; 1 = normal
species = 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
noiseType = 1;  % 1 for variable fGn (0 for fixed fGn)
implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% period histogram parameters

num_phbins = 32; % number of time bins in the the period histogram
phbinwidth = 1/F0/num_phbins; % binwidth in seconds;
t_st = 10e-3; % start time of analysis window for psth
ind_st = round(t_st*Fs); % index to start of analysis window for psth
t_end = floor((stimdur-20e-3)/T0)*T0+t_st; % end time of analysis window for psth; make it an integer multiple of the stimulus period
ind_end =round(t_end*Fs); % index to start of analysis window for psth

% generate stimulus waveform
pin(1:stimpts) = sqrt(2)*20e-6*sin(2*pi*F0*t); % unramped stimulus
pin(1:irpts) = pin(1:irpts).*(0:(irpts-1))/irpts; % apply onset ramp
pin(stimpts-irpts:stimpts)=pin(stimpts-irpts:stimpts).*(irpts:-1:0)/irpts; % apply offset ramp

stim_ph = zeros(numdbs,num_phbins); %initilaize array for stimulus waveform over one period
periodhistograms = zeros(numdbs,num_phbins); %initilaize array for period histogram
periodhistograms_scaled = zeros(numdbs,num_phbins); %initilaize array for period histogram

meanrates = zeros(1,numdbs); % create array for mean rate values at each stimulus level
SIs = zeros(1,numdbs); % create array for synchronization index values at each stimulus level
numspikes = zeros(1,numdbs); % create array for number of spikes in the period histogram at each stimulus level
vihcmaxs = zeros(1,numdbs); % create array for maximum vihc values at each stimulus level
phasediffs = zeros(1,numdbs); % create array for phase difference between period histogram and sine phase at each stimulus level
shiftvals = zeros(1,numdbs); % create array for shift in period histogram to bring it to sine phase at each stimulus level

h_io = figure; % create a figure object for the I/O function plots

for dblp = 1:length(stimdbs)

    disp(['Run ' num2str(dblp) '/' num2str(length(stimdbs))])

    stimdb = stimdbs(dblp); % stimulus level for this loop

    vihc = model_IHC_BEZ2018a(pin*10^(stimdb/20),CF,nrep_stim,1/Fs,sim_dur,cohc,cihc,species); % vihc output from model

    vihcmaxs(dblp) = max(abs(vihc)); % store maximum absolute value of vihc

    psth = model_Synapse_BEZ2018a(vihc,CF,nrep_stim,1/Fs,noiseType,implnt,spont,tabs,trel,expliketype); % psth output ready for raster plot

    meanrates(dblp) = sum(psth(ind_st:ind_end)) /(t_end - t_st)/ nrep_stim; % firing rates at each stim summed

    psth_ph_prep = psth(ind_st:ind_end); % psth from t_st to t_end, used to generated the period histogram
    tpsth_ph_prep = (0:(length(psth_ph_prep)-1))/Fs; % time array for psth being used to generated the period histogram
    
    for lp = 1:length(psth_ph_prep)
        phbin = round(rem(2*pi*F0*tpsth_ph_prep(lp),2*pi)/(2*pi*F0)/phbinwidth)+1; % find the period histogram bin in which the psth bin falls
        if phbin > num_phbins
            phbin = mod(phbin, num_phbins); % wrap around spike that occur after final period hisogram bin due to rounding
        end

        periodhistograms(dblp,phbin) = periodhistograms(dblp,phbin)+psth_ph_prep(lp); % put the psth bin value into the period histogram

    end

    SI_sin = periodhistograms(dblp,:) * sin(2*pi*(1:num_phbins)/num_phbins)'; % calculate the synchronization index in sine phase
    SI_cos= periodhistograms(dblp,:) * cos(2*pi*(1:num_phbins)/num_phbins)'; % calculate the synchronization index in cos phase

    SI_sin = SI_sin/sum(periodhistograms(dblp,:)); % normalize the synchronization index in sine phase
    SI_cos = SI_cos/sum(periodhistograms(dblp,:)); % normalize the synchronization index in cos phase

    SIs(dblp) = sqrt(SI_sin^2 + SI_cos^2); % calculate the synchronization index

    numspikes(dblp) = sum(periodhistograms(dblp,:)); % compute number of spikes for the period histogram at this stimulus level

    % [M_ph,ph_ph,f_ph] = fourier_dt(periodhistograms(dblp,:),1/phbinwidth,'half'); % calculate the phase of the period histogram
    ph_ph = unwrap(angle(fft(periodhistograms(dblp,:)))); % calculate the phase of the period histogram

    phasevals = (0:num_phbins-1)/num_phbins*(2*pi); % phase values in radians for each period histogram bin

    stim_ph(dblp,:) = sqrt(2)*1e3*20e-6*10^((stimdb)/20)*sin((phasevals)); % stimulus at sine phase for each period histogram bin in millipascal

    % [Mx_st,ph_st,f_st] = fourier_dt(stim_ph(dblp,:),1/phbinwidth,'half'); % calculate the phase of the stimulus
    ph_st = unwrap(angle(fft(stim_ph(dblp,:)))); % calculate the phase of the period histogram

    phasediffs(dblp) = ph_ph(2)-ph_st(2); % calc phase difference between the period histogram and the stimulus

    shiftvals(dblp) = round(phasediffs(dblp)/phasevals(2)); % compute how many time bins to shift the period histogram at this stimulus level

    periodhistograms_scaled(dblp,:) = periodhistograms(dblp,:)/nrep_stim/((t_end-t_st)/num_phbins); % scale the period histogram to be in spikes/s

    figure
    bar(phasevals,circshift(periodhistograms_scaled(dblp,:),shiftvals(dblp))) % psth in spikes/sec
    xlabel('Phase (rad)')
    ylabel('Instantaneous Rate (spikes/s)')
    title([explikestr '; ' spontstr '; ' num2str(stimdbs(dblp)) ' dB SPL; ' num2str(numspikes(dblp)) ' spikes'])
    hold on
    plot(phasevals,(max(periodhistograms_scaled(dblp,:))-min(periodhistograms_scaled(dblp,:)))/2*(1+sin(phasevals))+min(periodhistograms_scaled(dblp,:)));

    figure(h_io);
    scatter(stim_ph(dblp,:), circshift(periodhistograms_scaled(dblp,:),shiftvals(dblp)))
    hold on
    xlabel('Instantaneous Pressure (mPa)')
    ylabel('Instantaneous Rate (spikes/s)')
    title([explikestr '; ' spontstr])

end

%%

figure(h_io)
lgstr = num2str(stimdbs');
lgstr = [lgstr repmat(' dB SPL',length(stimdbs),1)];
legend(lgstr,'Location','northwest')

%%

figure
yyaxis right
plot(stimdbs,meanrates)
ylabel('Mean rate (spikes/s)')
yyaxis left
plot(stimdbs,SIs)
ylim([0 1])
ylabel('SI')
xlabel('Simulus Level (dB SPL)')
title([explikestr '; ' spontstr])


