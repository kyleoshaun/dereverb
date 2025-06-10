function [metrics, nalr_struct] = run_all_metrics(refstim, teststim, Fs_stim, HL, stimdb, h_ha_24k)
% Computes several monaural predictors of speech intelligibility (SI) and speech quality (SQ)
% Note: A NAL-R filter is applied to the test stimuli to compensate the impact
%       of hearing loss on audibility (i.e., to focus our attention on degraded
%       representations which are less easy to compensate)  
%
% Syntax:
%   [metrics, nalr_struct] = run_all_metrics(refstim, teststim, Fs_stim, HL, stimdb)
%
% Inputs:
% - refstim: Reference Signal in Pa (Clean, undistorted, unprocessed)
% - teststim: Test signal in Pa at the ear drum (Distorted, processed)
% - Fs_stim: Sample rate for stimuli (Hz)
% - HL: Hearing loss vector. dB loss at 6 audiometric frequencies: [250 500 1000 2000 4000 6000] Hz
% - stimdb: Acoustic level of clean stimulus (dB SPL). 
% - h_ha_24k: Hearing Aid gain filter (FIR) applied to the test signal for a sample rate of fs = 24 kHz
%   Note: Both test signals and the reference signals should be pre-calibrated 
% %       to the specified dB SPL.
%         If test signals include noise etc, the signal level should be 
%         calibrated before addition of these components
%
% Outputs:
% - metrics: All SI and SQ metrics 
%            SI: MR/FT NSIM, STMI, HASPI, STOI
%            SQ: HASQI, VISQOL
% - nalr_struct: Structure containing filter coefficients and other info about
%                the NAL-R filter (linear frequency selective gains, linear phase)
%                used to compensate hearing loss.

%% Compensate hearing loss: Apply NAL-R EQ with only linear gain compensation
% NOTES/QUESTIONS
% - This is usually run as part of HASQI (but not HASPI or any other metrics)
% - Done as preprocessing for all other metrics to neglect audibility


% Resample to 24 kHz (for auditory modeling in hearing aid gain (Matching fs used in NAL-R design)
[teststim_24k, fsamp] = eb_Resamp24kHz(teststim,  Fs_stim);

nsamp = length(teststim_24k);

% nfir = 140; %Length in samples of the FIR NAL-R EQ filter (24-kHz rate)
% [h_nalr,~]=eb_NALR(HL, nfir, fsamp); %Design the NAL-R filter

% Overwrite with custom hearing aid gain filter
nfir   = length(h_ha_24k) - 1;


% Apply the NAL-R filter
teststim_nalr_24k = conv(teststim_24k,  h_ha_24k);

% Apply Delay compensation
teststim_nalr_24k = teststim_nalr_24k(nfir+1:nfir+nsamp);

% Re-sample from 24kHz back to original stimulus sample rate
teststim_nalr = resample(teststim_nalr_24k,  Fs_stim, fsamp);


% Save NAL-R Filter Parameters
nalr_struct.h_nalr             = h_ha_24k;
nalr_struct.teststim_pre_nalr  = teststim;
nalr_struct.teststim_post_nalr = teststim_nalr;

%% TEMP

% Nfft = 4096;
% Nfft_welch  = 2^(nextpow2(length(teststim)));
% Nfft_welch  = min(Nfft_welch, Nfft);
% Nwin_welch  = Nfft_welch / 2;
% Nover_welch = Nfft_welch / 4;
% 
% [T_nalr, freqs_welch] = pwelch(teststim_nalr,  Nwin_welch, Nover_welch, Nfft_welch, Fs_stim);
% 
% [T,  ~] = pwelch(teststim,  Nwin_welch, Nover_welch, Nfft_welch, Fs_stim);
% 
% figure()
% subplot(2,1,1)
% semilogx(freqs_welch, 10*log10(T))
% xlabel('Frequency [Hz]')
% ylabel('dB')
% title('Test Stimulus before NAL-R Filter')
% subplot(2,1,2)
% semilogx(freqs_welch, 10*log10(T_nalr))
% xlabel('Frequency [Hz]')
% ylabel('dB')
% title('Test Stimulus after NAL-R Filter')
% 
% figure()
% freqz(h_nalr, 1, Nfft, Fs_stim);

%% Compute NSIM/STMI

[STMI, NSIM_FT, NSIM_MR] = nsim_and_stmi(teststim_nalr, refstim, Fs_stim, HL);

% Bypass to speed up execution for debug purposes
%STMI = 1;
%NSIM_FT = 1;
%NSIM_MR = 1;

%% Compute HASPI

[HASPI, ~]         = HASPI_v2(refstim, Fs_stim, teststim_nalr, Fs_stim, HL, stimdb);

%% Compute HASQI

[HASQI,~,~,~]         = HASQI_v2(refstim,Fs_stim,teststim_nalr, Fs_stim,HL,1,stimdb);

%% Compute STOI

STOI         = stoi(teststim, refstim, Fs_stim);
% Note: doesn't use NAL-R filtered test stimulus (No hearing loss included)


%% Compute VISQOL

VISQOL         = visqol(teststim, refstim, Fs_stim);
% Note: doesn't use NAL-R filtered test stimulus (No hearing loss included)


%% Save metrics to output

metrics.STMI           = STMI;
metrics.NSIM_FT        = NSIM_FT;
metrics.NSIM_MR        = NSIM_MR;
metrics.HASPI          = HASPI;
metrics.STOI           = STOI;
metrics.HASQI          = HASQI;
metrics.VISQOL         = VISQOL;

end