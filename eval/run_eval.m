function [metrics] = run_eval(source_data, noise_data, interf_data, fs, N60, s1_enable, s1_on_clean_speech, HL_NH, HL_HI, h_ha_NH, h_ha_HI, results_dir)
% Run DAP-dereverb algorithm for specified speech/noise/interference/reverb conditions
% and compute all SI/SQ metrics for the resulting unprocessed/processed signals.
%
% Syntax:
%   [metrics] = run_eval(source_data, noise_data, interf_data, fs, N60, s1_enable, s1_on_clean_speech, HL_NH, HL_HI, h_ha_NH, h_ha_HI, results_dir)
%
% Inputs:
% - source_data.signal:   Clean Speech Signal
% - source_data.rir_data: M-channel RIR data for source (cols = RIRs, M cols) 
% - source_data.memo:     Descriptor string (for file naming)
% - source_data.rir_memo: Short Descriptor for room (for file naming)
% - source_data.rir_desc: Long Descriptor for room (for log)
% - source_data.stimdb:   Source signal level to at the ear drum in dB SPL
%                         This function will automatically calibrate the 
%                         source signal level accordingly
% 
% - noise_data.enable:    Enable Noise
% - noise_data.signals:   M-channel noise signal data (cols = signals, M cols) -- cols length must match clean speech
% - noise_data.SNR_dB:       Signal-to-noise ratio (dB)
% - noise_data.memo:      Descriptor string (for file naming)
%
% - interf_data.enable:   Enable interference signal (secondary talker)
% - interf_data.signal:   Secondary clean speech signal (length must match clean speech)
% - interf_data.rir_data: M-channel RIR data for interfering talker 
% - interf_data.SIR_dB:      Signal-to-interference ratio (dB)
% - interf_data.memo:     Descriptor string (for file naming)
%
% - fs:  Sample rate for stimuli (Hz)
% - N60: T60 reverb time in samples (N60 = T60 * fs) for source_data.rir_data
% - s1_enable: 1 to enable the source-whitening stage in DAP-dereverb algo
% - s1_on_clean_speech: 1 to enable computation of DAP source whitening filter using clean speech (not blind)
% - HL_NH: Hearing loss vector for normal hearing listener (typically all zeros). dB loss at 6 audiometric frequencies: [250 500 1000 2000 4000 6000] Hz
% - HL_HI: Hearing loss vector for hearing impaired listener. dB loss at 6 audiometric frequencies: [250 500 1000 2000 4000 6000] Hz
% - h_ha_NH: Hearing Aid gain (FIR Filter, fs = 24 kHz) for NH listener
% - h_ha_HI: Hearing Aid gain (FIR Filter, fs = 24 kHz) for HI listener
% - run_memo: Descriptor for test/run which is added to plot titles
% - results_dir: Directory for saving results
%
% Outputs:
% - metrics.metrics_unproc_NH: All SI/SQ metrics for HL=HL_NH before DAP-dereverb (reverberant/unprocessed)
% - metrics.metrics_proc_NH:   All SI/SQ metrics for HL=HL_NH after DAP-dereverb (dereverberated/processed)
% - metrics.metrics_unproc_HI: All SI/SQ metrics for HL=HL_HI before DAP-dereverb (reverberant/unprocessed)
% - metrics.metrics_proc_HI:   All SI/SQ metrics for HL=HL_HI after DAP-dereverb (dereverberated/processed)
% - metrics.metrics_physical_unproc:   All physical metrics (not perceptual, e.g., DRR)
% - metrics.metrics_physical_proc:     All physical metrics (not perceptual, e.g., DRR)
%
% Note:
% - Channel Indexing convention: ch1 = Left1, ch2 = right1, ch3=left2, ch4=right2, ...

%%  Calling Script must include

% code_folder = "../../"; % Change this
% 
% addpath(sprintf('%s/dereverb/eval/', code_folder))
%
% addpath(sprintf('%s/dereverb/LPC/', code_folder))
% addpath(sprintf('%s/dereverb/utilities/matlab', code_folder))
% addpath(sprintf('%s/metrics/', code_folder))
% 
% addpath(sprintf('%s/samples', code_folder))
% addpath(sprintf('%s/RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/', code_folder))
% addpath(sprintf('%s/RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/', code_folder))
% addpath(sprintf('%s/RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/courtyard/', code_folder))
% 
% % Integration of EC with BEZ model
% addpath(sprintf("%s/metrics/BSIM2020_EC_BEZ2018a_model/", code_folder))
% 
% % HASPI/HASQI
% addpath(sprintf("%s/metrics/HASPIv2_HASQIv2_HAAQIv1_Common/", code_folder))
% 
% % BSIM
% addpath(sprintf('%s/metrics/BSIM_2020/', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/vad/', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/PreProc/', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/gammatonegram/', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/auditory/', code_folder))
% 
% addpath(genpath(sprintf('%s/metrics/BSIM_2020/Anechoic/', code_folder)))
% addpath(sprintf('%s/metrics/BSIM_2020/functions/', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/functions/gammatonefilter/', code_folder))
% addpath(sprintf('%s/metrics/BSIM_2020/functions/ltfat/', code_folder))
% 
% % Model (BEZ2018a)
% addpath(sprintf('%s/metrics/BEZ2018a_model/', code_folder))
% 
% % Other
% addpath(sprintf('%s/metrics/NSIM_Papers_ExampleCode 2/NSIMforBEZ2018a_model/', code_folder))
% addpath(sprintf('%s/metrics/STMI_Papers_ExampleCode 1/STMIforBEZ2018a_model', code_folder))

%% Parameters

% Extract Source Data
s              = source_data.signal;
G_source       = source_data.rir_data;
[L_channel, M] = size(G_source);
source_memo    = source_data.memo;
rir_memo       = source_data.rir_memo;
rir_desc       = source_data.rir_desc;

noise_memo = "NONE";
if noise_data.enable
    snr_db               = noise_data.SNR_dB;
    noise_memo           = sprintf("%s%.0fdB", noise_data.memo, snr_db);
else
    snr_db = Inf;
end

interf_memo = "NONE";
if interf_data.enable
    sir_db        = interf_data.SIR_dB;
    interf_memo   = sprintf("%s%.0fdB", interf_data.memo, sir_db);
else
    sir_db = Inf;
end

% Prediction Orders
T60_max = 500 / 1000;
N60_max = T60_max * fs;
p2 = round(N60_max / (M-1)); % Stage 2 MC-LPC order (Fixed)
%p2 = round(N60 / (M-1)); % Stage 2 MC-LPC order (Adaptive)
p1 = round(1.25*p2*(M-1)); % Stage 1 Source Whitening order

adjust_gain_ambiguity   = true;

% FFT/PSD params for plots
Nfft = 4096;
Nfft_welch  = 2^(nextpow2(length(s)));
Nfft_welch  = min(Nfft_welch, Nfft);
Nwin_welch  = Nfft_welch / 2;
Nover_welch = Nfft_welch / 4;

[Sm, w_welch] = pwelch(s, Nwin_welch, Nover_welch, Nfft_welch);
freqs_welch   = w_welch .* (fs / (2*pi));

%% Compute Microphone Signals

[Y, ~] = compute_mic_signals(source_data, noise_data, interf_data, fs);

%% Run Delay-and-Predict

T60_msec = (N60/fs) * 1000;

% Create Results Folder
s1_memo = "1s1";
if s1_enable == false
    s1_memo = "0s1";
end
parent_results_dir = results_dir;
results_dir = sprintf('%s/dap_rir%s_source%s_noise%s_interf%s_%.0fM_%s_%.0fT60_%.0fp1_%.0fp2_%.0fkFs_SNR%ddB_%s', parent_results_dir, rir_memo, source_memo, noise_memo, interf_memo, M, s1_memo, T60_msec, p1, p2, fs / 1000, snr_db, datetime('today'));
if ~exist(results_dir,'Dir')
    mkdir(results_dir)
end

% Start Logging Console if exist(
diary_fullpath = sprintf('%s/console.log', results_dir);
if exist(diary_fullpath, 'file')
    delete(diary_fullpath)
end
diary(diary_fullpath)

audiowrite(sprintf('%s/s.wav', results_dir), s .* (0.5 / max(abs(s))), fs);

fprintf("\nTest Conditions:\n")
fprintf(" - Source Signal = %s\n", source_data.memo)
fprintf(" - Source length = %d\n", length(source_data.signal))
fprintf(" - RIR = %s\n", source_data.rir_desc)
fprintf(" - RIR length = %.0f\n", L_channel)
fprintf(" - T60 = %.0f msec (N60 = %.0f samples)\n", T60_msec, N60)
fprintf(" - Noise Enabled = %d\n", noise_data.enable)
fprintf(" - SNR = %.0f dB\n", noise_data.SNR_dB)
fprintf(" - Noise Signal = %s\n", noise_data.memo)
fprintf(" - SIR = %.0f dB\n", interf_data.SIR_dB)
fprintf(" - Interference Signal = %s\n", interf_data.enable)

[s_est, dap_equalizer_struct] = delay_and_predict(Y, s, p1, p2, fs, s1_enable, s1_on_clean_speech, results_dir);

% Adjust for gain/scaling ambiguity
if adjust_gain_ambiguity
    s_est = (rms(s) / rms(s_est)) .* s_est;
end

H      = dap_equalizer_struct.H;
delays = dap_equalizer_struct.delays;
h0     = dap_equalizer_struct.h0;
delay_comp = dap_equalizer_struct.delay_comp;

% Save audio
audiowrite(sprintf('%s/s_training.wav',     results_dir), s .* (0.5 / max(abs(s))), fs);
audiowrite(sprintf('%s/y_training.wav',     results_dir), Y(:,1) .* (0.5 / max(abs(Y(:,1)))), fs);
audiowrite(sprintf('%s/s_est_training.wav', results_dir), s_est .* (0.5 / max(abs(s_est))), fs);

% [Tm, w_welch] = pwelch(s, Nwin_welch, Nover_welch, Nfft_welch);
% [Tm2, w_welch] = pwelch(Y(:,1), Nwin_welch, Nover_welch, Nfft_welch);
% [Tm3, w_welch] = pwelch(s_est, Nwin_welch, Nover_welch, Nfft_welch);
% figure()
% subplot(3,1,1)
% plot(10*log10(Tm))
% title("Clean Source")
% subplot(3,1,2)
% plot(10*log10(Tm2))
% title("Reverberant")
% subplot(3,1,3)
% plot(10*log10(Tm3))
% title("Dereverb")

%% MINT: Ideal solution based on known RIRs

run_MINT = 0;
if run_MINT
    fprintf("Running MINT (USING SAME L_g = p2+1) ... ");
    tau = 0;
    
    H_mint = MINT(G_source, tau, fs, p2+1, results_dir);
    
    save(sprintf('%s/H_mint.mat', results_dir),'H_mint');

    fprintf("Done\n\n")

end


%% Plots

%s_delayed  = filter([zeros(delay_comp,1) ; 1], 1, s);

%s_est = apply_dap_equalizer(dap_equalizer_struct, Y);
%s_est = (rms(s) / rms(s_est)) .* s_est;

e_est_dap = s - s_est;

Sm          = pwelch(s,         Nwin_welch, Nover_welch, Nfft_welch);
S_est_dap_m = pwelch(s_est,     Nwin_welch, Nover_welch, Nfft_welch);
E_est_dap_m = pwelch(e_est_dap, Nwin_welch, Nover_welch, Nfft_welch);

ylim_max = max([max(abs(s)), max(abs(s_est)), max(abs(e_est_dap))]) * 1.2;
figure()
subplot(3,1,1)
plot(s)
ylim([-1*ylim_max, ylim_max])
title('Source Signal')
subplot(3,1,2)
plot(s_est)
ylim([-1*ylim_max, ylim_max])
title('Estimated Source Signal (Delay-and-Predict)')
subplot(3,1,3)
plot(e_est_dap)
%ylim([-1*ylim_max, ylim_max])
title('Error')
sgtitle('Stage 2: MC-LPC Results')

saveas(gcf, sprintf('%s/S2_MCLPC_TD.fig', results_dir));
my_export_gcf(sprintf('%s/S2_MCLPC_TD.eps', results_dir), "TD_SIGNAL")

figure()
subplot(2,1,1)
plot(freqs_welch ./ 1000, 10*log10(Sm))
hold on;
plot(freqs_welch ./ 1000, 10*log10(S_est_dap_m))
legend('Source', 'Estimated Source')
title('Magnitude Spectrum')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(freqs_welch ./ 1000, 10*log10(E_est_dap_m))
title('Error Spectrum')
xlabel('Frequency [kHz]')
ylabel('dB')
sgtitle('Stage 2: MC-LPC Results')

saveas(gcf, sprintf('%s/S2_MCLPC_FD.fig', results_dir));
my_export_gcf(sprintf('%s/S2_MCLPC_FD.eps', results_dir), "FD_SIGNAL")

%% Generate/Plot IR for total DAP system (pass impulse through) and for DAP system alone

eir = apply_dap_equalizer(dap_equalizer_struct, G_source);

S_imp    = zeros(length(G_source(:,1)),M);
for ch = 1:M
    S_imp(1,ch) = 1;
end
ir_dap   = apply_dap_equalizer(dap_equalizer_struct, S_imp);

[G_channel_1, freqs_freqz] = freqz(G_source(:,1), 1, Nfft, fs);
Gm_channel_1 = abs(G_channel_1);
Gp_channel_1 = angle(G_channel_1);

G_channel_2  = freqz(G_source(:,2), 1, Nfft);
Gm_channel_2 = abs(G_channel_2);
Gp_channel_2 = angle(G_channel_2);

EIRft = fft(eir, Nfft);
EIRm  = abs(EIRft(1:(Nfft/2)));
EIRp  = angle(EIRft(1:(Nfft/2)));
freqs_fft = (0:(length(EIRm)-1))' .* (fs / Nfft);

IR_DAPft = fft(ir_dap, Nfft);
IR_DAPm  = abs(IR_DAPft(1:(Nfft/2)));
IR_DAPp  = angle(IR_DAPft(1:(Nfft/2)));
freqs_fft = (0:(length(IR_DAPm)-1))' .* (fs / Nfft);

ylim_min = min((max(eir(:,1))*-0.05), (min(eir(:,1))*1.1));
figure()
subplot(2,1,1)
plot((0:(length(G_source(:,1))-1)) .* (1/fs), G_source(:,1));
xlim([((length(G_source(:,1))/fs) * -0.02) Inf])
ylim([(min(G_source(:,1))*1.1) (max(G_source(:,1))*1.1)])
xlabel('Time [sec]')
title("Channel 1 Impulse Response")
subplot(2,1,2)
plot((0:(length(eir)-1)) .* (1/fs), eir);
xlim([((length(eir(:,1))/fs) * -0.02) Inf])
ylim([ylim_min (max(eir(:,1))*1.1)])
xlabel('Time [sec]')
title("Equalized Impulse Response")

saveas(gcf, sprintf('%s/EIR.fig', results_dir));
my_export_gcf(sprintf('%s/EIR.eps', results_dir), "TD_SIGNAL")


figure()
subplot(1,2,1)
plot(freqs_freqz ./ 1000, 20*log10(Gm_channel_1));
title('Channel 1 magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(1,2,2)
plot(freqs_freqz ./ 1000, Gp_channel_1);
title('Channel 1 phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

saveas(gcf, sprintf('%s/RTF.fig', results_dir));
my_export_gcf(sprintf('%s/RTF.eps', results_dir), "FD_SIGNAL")

figure()
subplot(2,1,1)
plot(freqs_fft ./ 1000, 20*log10(EIRm));
ylim([(20*log10(mean(EIRm)) - 3) (20*log10(mean(EIRm)) + 3)])
title('Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(freqs_fft ./ 1000, unwrap(EIRp));
ylim([-1 1])
title('Equalized phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

saveas(gcf, sprintf('%s/Equalized_RTF.fig', results_dir));
my_export_gcf(sprintf('%s/Equalized_RTF.eps', results_dir), "FD_SIGNAL")

figure()
subplot(2,1,1)
plot((0:(length(ir_dap)-1)) .* (1/fs), ir_dap);
xlim([((length(ir_dap)/fs) * -0.02) Inf])
ylim([(min(ir_dap)*1.1) (max(ir_dap)*1.1)])
xlabel("Full Delay-and-Predict Equalizer Impulse Response")
subplot(2,1,2)
plot(freqs_fft ./ 1000, 20*log10(IR_DAPm));
title('Full Delay-and-Predict Equalizer Magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')

saveas(gcf, sprintf('%s/dap_ir_mr.fig', results_dir));
my_export_gcf(sprintf('%s/dap_ir_mr.eps', results_dir), "FD_SIGNAL")

%% Energy Decay Curve Comparison

edc_rir_1 = EDC(G_source(:,1));
edc_dap   = EDC(eir);


figure()
plot((0:(length(G_source(:,1))-1)) .* (1/fs), 10*log10(edc_rir_1));
hold on;
plot((0:(length(G_source(:,1))-1)) .* (1/fs), 10*log10(edc_dap .* (max(edc_rir_1) / max(edc_dap))));
ylim([-70 6])
xlim([((length(G_source(:,1))/fs) * -0.02) Inf])
%xlim([0 (length(h_channel_1) * (1/fs))])
xlabel('Time [sec]')
ylabel('Energy remaining [dB]')
legend('Reverb Energy', 'Equalized Reverb Energy')
title('Energy Decay Curve')

saveas(gcf, sprintf('%s/EDC.fig', results_dir));
my_export_gcf(sprintf('%s/EDC.eps', results_dir), "TD_SIGNAL")

%% Spectrogram Results

% Spectrogram Parameters
N_window  = 256;
N_overlap = N_window / 2;
N_fft     = N_window;

% Compute Spectrograms
window = hamming(N_window);
[S_spec, f_S, t_S] = spectrogram(s, window, N_overlap, N_fft, fs);
[Y_spec, f_Y, t_Y] = spectrogram(Y(:,1), window, N_overlap, N_fft, fs);
[S_est_spec, f_S_est, t_S_est] = spectrogram(s_est, window, N_overlap, N_fft, fs);

figure()
subplot(3,1,1)
imagesc(t_S*1000, f_S/1e3, 20*log10(abs(S_spec)/sum(window)*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
xlabel('Time [msec]')
ylabel('Frequency [kHz]')
title('Clean Speech')
subplot(3,1,2)
imagesc(t_Y*1000, f_Y/1e3, 20*log10(abs(Y_spec)/sum(window)*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
xlabel('Time [msec]')
ylabel('Frequency [kHz]')
title('Reverberant Speech')
subplot(3,1,3)
imagesc(t_S_est*1000, f_S_est/1e3, 20*log10(abs(S_est_spec)/sum(window)*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
xlabel('Time [msec]')
ylabel('Frequency [kHz]')
title('De-reverberated Speech')
sgtitle('Dereverberation Training Results')

saveas(gcf, sprintf('%s/Spectrogram.fig', results_dir));
my_export_gcf(sprintf('%s/Spectrogram.eps', results_dir), "SPECTROGRAM")

%% Extra plots

enable_extra_plots = 1;
if enable_extra_plots

% NOTE: Assumes M = 4
if M ~= 4
    error("Extra plots only work for M=4")
end

g1 = G_source(:,1);
g2 = G_source(:,2);
g3 = G_source(:,3);
g4 = G_source(:,4);

g1 = filter([zeros(delays(1),1) ; 1], 1, g1);
g2 = filter([zeros(delays(2),1) ; 1], 1, g2);
g3 = filter([zeros(delays(3),1) ; 1], 1, g3);
g4 = filter([zeros(delays(4),1) ; 1], 1, g4);

h_pred_1_from_1 = H(:,1);
h_pred_1_from_2 = H(:,2);
h_pred_1_from_3 = H(:,3);
h_pred_1_from_4 = H(:,4);

h_pred_2_from_1 = H(:,5);
h_pred_2_from_2 = H(:,6);
h_pred_2_from_3 = H(:,7);
h_pred_2_from_4 = H(:,8);

h_pred_3_from_1 = H(:,9);
h_pred_3_from_2 = H(:,10);
h_pred_3_from_3 = H(:,11);
h_pred_3_from_4 = H(:,12);

h_pred_4_from_1 = H(:,13);
h_pred_4_from_2 = H(:,14);
h_pred_4_from_3 = H(:,15);
h_pred_4_from_4 = H(:,16);

eir_1 = g1 - ((filter(h_pred_1_from_1, 1, g1)) + ...
              (filter(h_pred_1_from_2, 1, g2)) + ...
              (filter(h_pred_1_from_3, 1, g3)) + ...
              (filter(h_pred_1_from_4, 1, g4)));

eir_2 = g2 - ((filter(h_pred_2_from_1, 1, g1)) + ...
              (filter(h_pred_2_from_2, 1, g2)) + ...
              (filter(h_pred_2_from_3, 1, g3)) + ...
              (filter(h_pred_2_from_4, 1, g4)));

eir_3 = g3 - ((filter(h_pred_3_from_1, 1, g1)) + ...
              (filter(h_pred_3_from_2, 1, g2)) + ...
              (filter(h_pred_3_from_3, 1, g3)) + ...
              (filter(h_pred_3_from_4, 1, g4)));

eir_4 = g4 - ((filter(h_pred_4_from_1, 1, g1)) + ...
              (filter(h_pred_4_from_2, 1, g2)) + ...
              (filter(h_pred_4_from_3, 1, g3)) + ...
              (filter(h_pred_4_from_4, 1, g4)));


% eir = eir_1 .* h0(1) + ...
%       eir_2 .* h0(2) + ...
%       eir_3 .* h0(3) + ...
%       eir_4 .* h0(4);

EIR_1  = freqz(eir_1, 1, Nfft);
EIR_2  = freqz(eir_2, 1, Nfft);
EIR_3  = freqz(eir_3, 1, Nfft);
EIR_4  = freqz(eir_4, 1, Nfft);

EIR_1_m = abs(EIR_1);
EIR_2_m = abs(EIR_2);
EIR_3_m = abs(EIR_3);
EIR_4_m = abs(EIR_4);

rir_ylim_max = max([max(g1) ; max(g2) ; max(g3) ; max(g4)]);
rir_ylim_min = min(1.1*min([min(g1) ; min(g2) ; min(g3) ; min(g4)]), rir_ylim_max*-0.05);

eir_ylim_max = max([max(eir_1) ; max(eir_2) ; max(eir_3) ; max(eir_4)]);
eir_ylim_min = min(1.1*min([min(eir_1) ; min(eir_2) ; min(eir_3) ; min(eir_4)]), eir_ylim_max*-0.05);

fig_size = [10 10 1300 600];

figure()
set(gcf,'Position',fig_size)

subplot(4,2,1)
plot((0:(length(g1)-1)) .* (1/fs), g1);
ylim([rir_ylim_min rir_ylim_max])
xlim([((length(g1)/fs) * -0.02) Inf])
xlabel('Time [sec]')
title("Impulse Response (Channel 1)")
subplot(4,2,2)
plot((0:(length(eir_1)-1)) .* (1/fs), eir_1);
xlim([((length(eir_1)/fs) * -0.02) Inf])
ylim([eir_ylim_min eir_ylim_max])
xlabel('Time [sec]')
title("Equalized Impulse Response (MC-LP Prediction of Channel 1)")

subplot(4,2,3)
plot((0:(length(g2)-1)) .* (1/fs), g2);
xlim([((length(g2)/fs) * -0.02) Inf])
ylim([rir_ylim_min rir_ylim_max])
xlabel('Time [sec]')
title("Impulse Response (Channel 2)")
subplot(4,2,4)
plot((0:(length(eir_2)-1)) .* (1/fs), eir_2);
xlim([((length(eir_2)/fs) * -0.02) Inf])
ylim([eir_ylim_min eir_ylim_max])
xlabel('Time [sec]')
title("Equalized Impulse Response (MC-LP Prediction of Channel 2)")

subplot(4,2,5)
plot((0:(length(g3)-1)) .* (1/fs), g3);
xlim([((length(g3)/fs) * -0.02) Inf])
ylim([rir_ylim_min rir_ylim_max])
xlabel('Time [sec]')
title("Impulse Response (Channel 3)")
subplot(4,2,6)
plot((0:(length(eir_3)-1)) .* (1/fs), eir_3);
xlim([((length(eir_3)/fs) * -0.02) Inf])
ylim([eir_ylim_min eir_ylim_max])
xlabel('Time [sec]')
title("Equalized Impulse Response (MC-LP Prediction of Channel 3)")

subplot(4,2,7)
plot((0:(length(g4)-1)) .* (1/fs), g4);
xlim([((length(g4)/fs) * -0.02) Inf])
ylim([rir_ylim_min rir_ylim_max])
xlabel('Time [sec]')
title("Impulse Response (Channel 1)")
subplot(4,2,8)
plot((0:(length(eir_4)-1)) .* (1/fs), eir_4);
xlim([((length(eir_4)/fs) * -0.02) Inf])
ylim([eir_ylim_min eir_ylim_max])
xlabel('Time [sec]')
title("Equalized Impulse Response (MC-LP Prediction of Channel 4)")

saveas(gcf, sprintf('%s/MC_EIR.fig', results_dir));
my_export_gcf(sprintf('%s/MC_EIR.eps', results_dir), "TD_SIGNAL")

figure()
subplot(2,2,1)
plot(freqs_freqz ./ 1000, 20*log10(EIR_1_m))
title('Equalized RTF (MC-LP of Channel 1)')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,2)
plot(freqs_freqz ./ 1000, 20*log10(EIR_2_m))
title('Equalized RTF (MC-LP of Channel 2)')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,3)
plot(freqs_freqz ./ 1000, 20*log10(EIR_3_m))
title('Equalized RTF (MC-LP of Channel 3)')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,4)
plot(freqs_freqz ./ 1000, 20*log10(EIR_4_m))
title('Equalized RTF (MC-LP of Channel 4)')
xlabel('Frequency [kHz]')
ylabel('dB')

saveas(gcf, sprintf('%s/MC_Equalized_RTF.fig', results_dir));
my_export_gcf(sprintf('%s/MC_Equalized_RTF.eps', results_dir), "FD_SIGNAL")

end % Extra plots

%% Test: Apply dereverb to a short speech signal SA2.wav
%  NOTE: A DAP-EQ that over-whitens SA1.WAV due to p1 being too low will
%        not over-whiten SA2.WAV, it will actually add a reverb-like effect
%        to SA2.WAV

[s,fs] = audioread("SA2.WAV");
source_data.signal = s;

% Run Metrics on reverb only
noise_data_noNoise.enable    = false;
noise_data_noNoise.signals   = [];
noise_data_noNoise.SNR_dB    = Inf;
noise_data_noNoise.memo      = "no noise";

interf_data_noInterf.enable   = false;
interf_data_noInterf.signal   = [0];
interf_data_noInterf.rir_data = [0];
interf_data_noInterf.SIR_dB   = Inf;

[Y, refstim] = compute_mic_signals(source_data, noise_data_noNoise, interf_data_noInterf, fs);

s_est = apply_dap_equalizer(dap_equalizer_struct, Y);

% Adjust for gain/scaling ambiguity
if adjust_gain_ambiguity
    s_est = (rms(refstim) / rms(s_est)) .* s_est;
end

% Save audio
audiowrite(sprintf('%s/s_test_SA2.wav',     results_dir), s .* (0.5 / max(abs(s))),      fs);
audiowrite(sprintf('%s/y_test_SA2.wav',     results_dir), Y(:,1) .* (0.5 / max(abs(Y(:,1)))), fs);
audiowrite(sprintf('%s/s_test_est_SA2.wav', results_dir), s_est .* (0.5 / max(abs(s_est))),  fs);

% Compute MINT Signal Results
% % s_est_MINT = apply_dap_equalizer(dap_equalizer_struct, Y);
% s_est_MINT = zeros(size(Y(:, 1)));
% for ch_idx = 1:M
%     y_ch       = Y(:, ch_idx);
%     h_mint_ch  = H_mint(:, ch_idx);
%     s_est_MINT = s_est_MINT + filter(h_mint_ch, 1, y_ch);
% end


% Spectrogram Parameters
N_window  = 256;
N_overlap = N_window / 2;
N_fft     = N_window;

% Compute Spectrograms
window = hamming(N_window);
[S_spec, f_S, t_S] = spectrogram(refstim, window, N_overlap, N_fft, fs);
[Y_spec, f_Y, t_Y] = spectrogram(Y(:,1), window, N_overlap, N_fft, fs);
[S_est_spec, f_S_est, t_S_est] = spectrogram(s_est, window, N_overlap, N_fft, fs);
%[S_est_MINT_spec, f_S_est_MINT, t_S_est_MINT] = spectrogram(s_est_MINT, window, N_overlap, N_fft, fs);

figure()
subplot(3,1,1)
imagesc(t_S*1000, f_S/1e3, 20*log10(abs(S_spec)/sum(window)*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
xlabel('Time [msec]')
ylabel('Frequency [kHz]')
title('Clean Speech')
subplot(3,1,2)
imagesc(t_Y*1000, f_Y/1e3, 20*log10(abs(Y_spec)/sum(window)*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
xlabel('Time [msec]')
ylabel('Frequency [kHz]')
title('Reverberant Speech')
subplot(3,1,3)
imagesc(t_S_est*1000, f_S_est/1e3, 20*log10(abs(S_est_spec)/sum(window)*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
xlabel('Time [msec]')
ylabel('Frequency [kHz]')
title('De-reverberated Speech')
sgtitle('Dereverberation Results (Reapplied DAP-EQ to SA2.wav, without noise)')

saveas(gcf, sprintf('%s/Spectrogram_reappliedToSA2.fig', results_dir));
my_export_gcf(sprintf('%s/Spectrogram_reappliedToSA2.eps', results_dir), "SPECTROGRAM")

% MINT Spectrogram Results
% figure()
% subplot(3,1,1)
% imagesc(t_S*1000, f_S/1e3, 20*log10(abs(S_spec)/sum(window)*sqrt(2)/20e-6));
% axis xy; axis tight;
% hcb = colorbar;
% set(get(hcb,'ylabel'),'string','SPL')
% xlabel('Time [msec]')
% ylabel('Frequency [kHz]')
% title('Clean Speech')
% subplot(3,1,2)
% imagesc(t_Y*1000, f_Y/1e3, 20*log10(abs(Y_spec)/sum(window)*sqrt(2)/20e-6));
% axis xy; axis tight;
% hcb = colorbar;
% set(get(hcb,'ylabel'),'string','SPL')
% xlabel('Time [msec]')
% ylabel('Frequency [kHz]')
% title('Reverberant Speech')
% subplot(3,1,3)
% imagesc(t_S_est_MINT*1000, f_S_est_MINT/1e3, 20*log10(abs(S_est_MINT_spec)/sum(window)*sqrt(2)/20e-6));
% axis xy; axis tight;
% hcb = colorbar;
% set(get(hcb,'ylabel'),'string','SPL')
% xlabel('Time [msec]')
% ylabel('Frequency [kHz]')
% title('De-reverberated Speech (MINT)')
% sgtitle('MINT Dereverberation Results (Reapplied DAP-EQ to SA2.wav, without noise)')
% 
% saveas(gcf, sprintf('%s/Spectrogram_MINT_reappliedToSA2.fig', results_dir));
% my_export_gcf(sprintf('%s/Spectrogram_MINT_reappliedToSA2.eps', results_dir), "SPECTROGRAM")


%% Compute SI/SQ Metrics (Without noise)

stimdb          = source_data.stimdb; % dB SPL
Fs_stim         = fs;
test_reverb     = Y(:,1);
test_processed  = s_est;

[metrics_unproc_NH, ~] = run_all_metrics(refstim, test_reverb,    Fs_stim, HL_NH, stimdb, h_ha_NH, results_dir);
[metrics_proc_NH,   ~] = run_all_metrics(refstim, test_processed, Fs_stim, HL_NH, stimdb, h_ha_NH, results_dir);

[metrics_unproc_HI, ~] = run_all_metrics(refstim, test_reverb,    Fs_stim, HL_HI, stimdb, h_ha_HI, results_dir);
[metrics_proc_HI,   ~] = run_all_metrics(refstim, test_processed, Fs_stim, HL_HI, stimdb, h_ha_HI, results_dir);

% audiowrite(sprintf('%s/ec_out_unproc.wav', results_dir), ec_out_unproc .* (0.5 / max(abs(ec_out_unproc))), fs);
% audiowrite(sprintf('%s/ec_out_proc.wav',   results_dir), ec_out_proc .* (0.5 / max(abs(ec_out_proc))),   fs);
% 
% audiowrite(sprintf('%s/test_stim_unproc_postNALR.wav', results_dir), ...
%            nalr_struct_unproc_HI.teststim_left_post_nalr .* ...
%            (0.5 / max(abs(nalr_struct_unproc_HI.teststim_left_post_nalr))), fs);
% audiowrite(sprintf('%s/test_stim_proc_postNALR.wav',   results_dir), ...
%            nalr_struct_proc_HI.teststim_left_post_nalr .* ...
%            (0.5 / max(abs(nalr_struct_proc_HI.teststim_left_post_nalr))), fs);

%% Physical metrics of reverb

C50_unproc = clarity(G_source(:,1), fs);
C50_proc   = clarity(eir, fs);

metrics_physical_unproc.C50 = C50_unproc;
metrics_physical_proc.C50   = C50_proc;

%% LOG RESULTS

fprintf("\n")
fprintf("Results for NH Listener:")
fprintf("\nunprocessed metric --> processed metric\n")
fprintf("SI RESULTS: -----------------------------------\n")
fprintf("HASPI         = %.8f --> %.8f\n", metrics_unproc_NH.HASPI,         metrics_proc_NH.HASPI);
fprintf("STOI          = %.8f --> %.8f\n", metrics_unproc_NH.STOI,          metrics_proc_NH.STOI);
fprintf("NSIM_FT       = %.8f --> %.8f\n", metrics_unproc_NH.NSIM_FT,       metrics_proc_NH.NSIM_FT);
fprintf("NSIM_MR       = %.8f --> %.8f\n", metrics_unproc_NH.NSIM_MR,       metrics_proc_NH.NSIM_MR);
fprintf("STMI          = %.8f --> %.8f\n", metrics_unproc_NH.STMI,          metrics_proc_NH.STMI);
fprintf("-----------------------------------------------")
fprintf("\n")
fprintf("SQ RESULTS: -----------------------------------\n")
fprintf("HASQI         = %.8f --> %.8f\n", metrics_unproc_NH.HASQI,          metrics_proc_NH.HASQI);
fprintf("VISQOL        = %.8f --> %.8f\n", metrics_unproc_NH.VISQOL,         metrics_proc_NH.VISQOL);
fprintf("-----------------------------------------------\n")
fprintf("\n")
fprintf("\n")
fprintf("Results for HI Listener:")
fprintf("\nunprocessed metric --> processed metric\n")
fprintf("SI RESULTS: -----------------------------------\n")
fprintf("HASPI         = %.8f --> %.8f\n", metrics_unproc_HI.HASPI,         metrics_proc_HI.HASPI);
fprintf("STOI          = %.8f --> %.8f\n", metrics_unproc_HI.STOI,          metrics_proc_HI.STOI);
fprintf("NSIM_FT       = %.8f --> %.8f\n", metrics_unproc_HI.NSIM_FT,       metrics_proc_HI.NSIM_FT);
fprintf("NSIM_MR       = %.8f --> %.8f\n", metrics_unproc_HI.NSIM_MR,       metrics_proc_HI.NSIM_MR);
fprintf("STMI          = %.8f --> %.8f\n", metrics_unproc_HI.STMI,          metrics_proc_HI.STMI);
fprintf("-----------------------------------------------")
fprintf("\n")
fprintf("SQ RESULTS: -----------------------------------\n")
fprintf("HASQI         = %.8f --> %.8f\n", metrics_unproc_HI.HASQI,          metrics_proc_HI.HASQI);
fprintf("VISQOL        = %.8f --> %.8f\n", metrics_unproc_HI.VISQOL,         metrics_proc_HI.VISQOL);
fprintf("-----------------------------------------------\n")
fprintf("\n")
fprintf("Physical Metrics Results ----------------------\n")
fprintf("\nunprocessed metric --> processed metric\n")
fprintf("Clarity (C50) = %.8f dB --> %.8f dB\n", 10*log10(metrics_physical_unproc.C50), 10*log10(metrics_physical_proc.C50));
fprintf("-----------------------------------------------\n")

%% Save results

metrics.metrics_unproc_NH = metrics_unproc_NH;
metrics.metrics_proc_NH   = metrics_proc_NH;

metrics.metrics_unproc_HI = metrics_unproc_HI;
metrics.metrics_proc_HI   = metrics_proc_HI;

metrics.metrics_physical_unproc = metrics_physical_unproc;
metrics.metrics_physical_proc   = metrics_physical_proc;

% Stop logging console
diary off;

end