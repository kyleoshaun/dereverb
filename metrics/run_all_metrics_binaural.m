function [metrics] = run_all_metrics_binaural(refstim, teststim_left, teststim_right, Fs_stim, HL, stimdb, h_ha_24k)
% Computes several monaural predictors of speech intelligibility (SI) and speech quality (SQ) with a binaural EC front-end
% Note: A NAL-R filter is applied to the test stimuli to compensate the impact
%       of hearing loss on audibility (i.e., to focus our attention on degraded
%       representations which are less easy to compensate)  
%
% Syntax:
%   [metrics, nalr_struct] = run_all_metrics(refstim, teststim, Fs_stim, HL, stimdb)
%
% Inputs:
% - refstim: Reference Signal in Pa (Clean, undistorted, unprocessed)
% - teststim_left: Left-ear est signal in Pa at the ear drum (Distorted, processed)
% - teststim_left: Right-ear est signal in Pa at the ear drum (Distorted, processed)
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

%% EC Front End

ec_out  = EC_FrontEnd(teststim_left, teststim_right, Fs_stim);  

% NOTE: Acoustic levels are not recalibrated after the EC algorithm, 
%       There is an assumption here that the gain applied by EC on the desired signal
%       is approximately unity. 

% Match ref signal to length
refstim             = [refstim ; 0];
teststim_left       = [teststim_left ; 0];
teststim_right      = [teststim_right ; 0];

%% Compensate hearing loss: Apply Linear hearing aid gain

%Resample to 24 kHz (for auditory modeling in hearing aid gain (Matching fs used in NAL-R design)
[ec_out_24k, fsamp] = eb_Resamp24kHz(ec_out,  Fs_stim);

nsamp = length(ec_out_24k);

% Overwrite with custom hearing aid gain filter
nfir   = length(h_ha_24k) - 1;

% Apply the hearing aid gain filter
ec_out_postGain_24k = conv(ec_out_24k,  h_ha_24k);

% Apply Delay compensation
ec_out_postGain_24k = ec_out_postGain_24k(nfir+1:nfir+nsamp);

% Re-sample from 24kHz back to original stimulus sample rate
ec_out_postGain = resample(ec_out_postGain_24k,  Fs_stim, fsamp);

% ec_out_24k          = resample(ec_out, 24000, Fs_stim);
% ec_out_postGain_24k = filter(h_ha_24k, 1, ec_out_24k);
% ec_out_postGain     = resample(ec_out_postGain_24k, Fs_stim, 24000);

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

%[STMI, NSIM_FT, NSIM_MR] = nsim_and_stmi(teststim_postGain, refstim, Fs_stim, HL);

STMI = 1;
NSIM_FT = 1;
NSIM_MR = 1;

%% Compute HASPI

[HASPI, ~]         = HASPI_v2(refstim, Fs_stim, ec_out_postGain, Fs_stim, HL, stimdb);

%% Compute HASQI

[HASQI,~,~,~]         = HASQI_v2(refstim,Fs_stim,ec_out, Fs_stim,HL,1,stimdb);
% Note: doesn't use NAL-R filtered test stimulus (already part of HASQI internally)


%% Compute STOI

STOI         = stoi(ec_out, refstim, Fs_stim);
% Note: doesn't use NAL-R filtered test stimulus (No hearing loss included)


%% Compute VISQOL

VISQOL         = visqol(ec_out, refstim, Fs_stim);
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