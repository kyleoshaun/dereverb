addpath('../../samples/')
addpath('../../AIR_Databases/AIR_1_4_BinauralDatabase/')
addpath('../../AIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/')
addpath('../../AIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/')

close all
clear
clc

% Shuffle and save seed for reproducibility
rng('shuffle')
seed = abs(round(randn(1,1) * 10000));
rng(seed);
rng(724); % Reproduce specific seed

%% Parameters

% Script Modes (Ordered by decreasing presedence)
test_stage2_alone       = false;
stage_1_on_clean_speech = false;
stage_1_on_8_mics       = false;
stage_1_on_4_mics       = false;

adjust_gain_ambiguity   = true;
apply_time_alignment    = false;
SNR_dB                  = 300;

% Choose Source signal
%[s_frame,fs] = audioread("assgn1.wav");
%[s_frame,fs] = audioread("SA1_phoneme.wav");
[s_frame,fs] = audioread("SA1.wav");

% Loop to synthetically increase length
% NOTE: MC-LPC stage performance is highly influenced by length of signal
%       Even if stage 1 is bypassed (For very long signals stage 1 bypassed
%       Results in nearly zero forcing)
num_loops_s_frame = round(8.8*fs / length(s_frame));
s_frame = repmat(s_frame, num_loops_s_frame, 1);

if test_stage2_alone
    % Set s_frame to white noise (rather than speech) to skip source whitening and
    % debug MC-LPC alone
    s_frame = randn(1000000,1); 
    % Note: The longer the channel IR, the more data is needed,
    % Length 5000000 will give near zero forcing for L_channel = 1000 (This is 312 sec at fs=16kHz though)
end
%s_frame = randn(500000,1); fs = 16000;

% FFT/PSD params for plots
Nfft = 4096;
Nfft_welch  = 2^(nextpow2(length(s_frame)));
Nfft_welch  = min(Nfft_welch, Nfft);
Nwin_welch  = Nfft_welch / 2;
Nover_welch = Nfft_welch / 4;

[Sm_frame, w_welch] = pwelch(s_frame, Nwin_welch, Nover_welch, Nfft_welch);
freqs_welch         = w_welch .* (fs / (2*pi));


%% Channel

% Exponential decay curve
L_channel = 4000;
tau = L_channel/4;
exp_decay = exp(-1 .* (1:L_channel)' ./ tau);

% Option 1: Generate random AIR
b_air_2 = randn(L_channel,1);
b_air_1 = randn(L_channel,1);
b_air_3 = randn(L_channel,1);
b_air_4 = randn(L_channel,1);
b_air_5 = randn(L_channel,1);
b_air_6 = randn(L_channel,1);
b_air_7 = randn(L_channel,1);
b_air_8 = randn(L_channel,1);

a_air_1 = 1;
a_air_2 = 1;
a_air_3 = 1;
a_air_4 = 1;
a_air_5 = 1;
a_air_6 = 1;
a_air_7 = 1;
a_air_8 = 1;

% Apply synthetic exponential decay to AIRs
b_air_1 = b_air_1(1:L_channel) .* exp_decay;
b_air_2 = b_air_2(1:L_channel) .* exp_decay;
if stage_1_on_4_mics || stage_1_on_8_mics
    b_air_3 = b_air_3(1:L_channel) .* exp_decay;
    b_air_4 = b_air_4(1:L_channel) .* exp_decay;
    b_air_5 = b_air_5(1:L_channel) .* exp_decay;
    b_air_6 = b_air_6(1:L_channel) .* exp_decay;
    b_air_7 = b_air_7(1:L_channel) .* exp_decay;
    b_air_8 = b_air_8(1:L_channel) .* exp_decay;
end

% Synthetic 2-sided AIR
%b_air_1 = [flip(b_air_1) ; b_air_1];
%b_air_2 = [flip(b_air_2) ; b_air_2];

% % Add Synthetic time delay
% delay_1 = 0;
% delay_2 = 0;
% b_air_1 = [zeros(delay_1,1) ; b_air_1 ; zeros(delay_2,1)];
% b_air_2 = [zeros(delay_2,1) ; b_air_2  ; zeros(delay_1,1)];

% Option 2: Real AIR
% Channels:
% - 1: Left Front
% - 2: Right Front
% - 3: Left Middle
% - 4: Right Middle
% - 5: Left Rear
% - 6: Front Rear
head_orientation = 2;
speaker_loc      = 'A';
data_set         = 'bte';
HRIR_data = loadHRIR('office_II', head_orientation, speaker_loc, data_set);
if HRIR_data.fs ~= fs
    [resample_p, resample_q] = rat(fs / HRIR_data.fs);
    HRIR_data.data = resample(HRIR_data.data, resample_p, resample_q);
end

b_air_1 = HRIR_data.data(:,1);
b_air_2 = HRIR_data.data(:,2);
b_air_3 = HRIR_data.data(:,3);
b_air_4 = HRIR_data.data(:,4);
b_air_5 = HRIR_data.data(:,5);
b_air_6 = HRIR_data.data(:,6);
b_air_7 = b_air_5;
b_air_8 = b_air_6;

trunc_length = 4000;
tof_bulk = 1;
b_air_1 = b_air_1(tof_bulk:trunc_length);
b_air_2 = b_air_2(tof_bulk:trunc_length);
b_air_3 = b_air_3(tof_bulk:trunc_length);
b_air_4 = b_air_4(tof_bulk:trunc_length);
b_air_5 = b_air_5(tof_bulk:trunc_length);
b_air_6 = b_air_6(tof_bulk:trunc_length);

tof_1 = 107;
tof_2 = 105;
tof_3 = 108;
tof_4 = 105;
tof_5 = 108;
tof_6 = 106;

% Zero out onset
% b_air_1(1:tof_1) = zeros(tof_1,1);
% b_air_2(1:tof_2) = zeros(tof_2,1);
% b_air_3(1:tof_3) = zeros(tof_3,1);
% b_air_4(1:tof_4) = zeros(tof_4,1);
% b_air_5(1:tof_5) = zeros(tof_5,1);
% b_air_6(1:tof_6) = zeros(tof_6,1);

% Manual time alignment via removal of delay
min_tof = min([tof_1, tof_2, tof_3, tof_4, tof_5, tof_6]);
delay_1 = tof_1 - min_tof;
delay_2 = tof_2 - min_tof;
delay_3 = tof_3 - min_tof;
delay_4 = tof_4 - min_tof;
delay_5 = tof_5 - min_tof;
delay_6 = tof_6 - min_tof;
b_air_1 = [b_air_1((delay_1+1):end) ; zeros(delay_1,1)];
b_air_2 = [b_air_2((delay_2+1):end) ; zeros(delay_2,1)];
b_air_3 = [b_air_3((delay_3+1):end) ; zeros(delay_3,1)];
b_air_4 = [b_air_4((delay_4+1):end) ; zeros(delay_4,1)];
b_air_5 = [b_air_5((delay_5+1):end) ; zeros(delay_5,1)];
b_air_6 = [b_air_6((delay_6+1):end) ; zeros(delay_6,1)];

L_channel = length(b_air_1);
p2 = round(L_channel * 1.25); % Stage 2 MC-LPC order
p1 = round(p2 + L_channel * 1.25); % Stage 1 Source Whitening order
%p1 = 200;

%% Compute Microphone Signals

% Compute reverberant signals
y1_frame = filter(b_air_1, a_air_1, s_frame);
y2_frame = filter(b_air_2, a_air_2, s_frame);
if stage_1_on_4_mics || stage_1_on_8_mics
    y3_frame = filter(b_air_3, a_air_3, s_frame);
    y4_frame = filter(b_air_4, a_air_4, s_frame);
    y5_frame = filter(b_air_5, a_air_5, s_frame);
    y6_frame = filter(b_air_6, a_air_6, s_frame);
    y7_frame = filter(b_air_7, a_air_7, s_frame);
    y8_frame = filter(b_air_8, a_air_8, s_frame);
end

% Add measurement noise
% NOTE: This is essentially regularizing (biasing) R_mc so dont make it unrealistically high
noise_mag = rms(s_frame) .* 10 ^ (-1*SNR_dB / 20);
y1_frame = y1_frame + noise_mag .* randn(length(y1_frame), 1);
y2_frame = y2_frame + noise_mag .* randn(length(y2_frame), 1);
if stage_1_on_4_mics || stage_1_on_8_mics
    y3_frame = y3_frame + noise_mag .* randn(length(y3_frame), 1);
    y4_frame = y4_frame + noise_mag .* randn(length(y4_frame), 1);
    y5_frame = y5_frame + noise_mag .* randn(length(y5_frame), 1);
    y6_frame = y6_frame + noise_mag .* randn(length(y6_frame), 1);
    y7_frame = y7_frame + noise_mag .* randn(length(y7_frame), 1);
    y8_frame = y8_frame + noise_mag .* randn(length(y8_frame), 1);
end

%% Interfering talker

[sv_frame,fs] = audioread("SA2.wav");

num_loops_s_frame = round(10*fs / length(sv_frame));
sv_frame = repmat(sv_frame, num_loops_s_frame, 1);

% Option 2: Real AIR
% Channels:
% - 1: Left Front
% - 2: Right Front
% - 3: Left Middle
% - 4: Right Middle
% - 5: Left Rear
% - 6: Front Rear
head_orientation = 2;
speaker_loc      = 'C';
data_set         = 'bte';
HRIR_data = loadHRIR('office_II', head_orientation, speaker_loc, data_set);
if HRIR_data.fs ~= fs
    [resample_p, resample_q] = rat(fs / HRIR_data.fs);
    HRIR_data.data = resample(HRIR_data.data, resample_p, resample_q);
end

b_air_v_1 = HRIR_data.data(:,1);
b_air_v_2 = HRIR_data.data(:,2);
b_air_v_3 = HRIR_data.data(:,3);
b_air_v_4 = HRIR_data.data(:,4);
b_air_v_5 = HRIR_data.data(:,5);
b_air_v_6 = HRIR_data.data(:,6);
b_air_v_7 = b_air_v_5;
b_air_v_8 = b_air_v_6;

b_air_v_1 = b_air_v_1(tof_bulk:trunc_length);
b_air_v_2 = b_air_v_2(tof_bulk:trunc_length);
b_air_v_3 = b_air_v_3(tof_bulk:trunc_length);
b_air_v_4 = b_air_v_4(tof_bulk:trunc_length);
b_air_v_5 = b_air_v_5(tof_bulk:trunc_length);
b_air_v_6 = b_air_v_6(tof_bulk:trunc_length);
b_air_v_7 = b_air_v_7(tof_bulk:trunc_length);
b_air_v_8 = b_air_v_8(tof_bulk:trunc_length);

% Compute reverberant signals
v1_frame = filter(b_air_v_1, 1, sv_frame);
v2_frame = filter(b_air_v_2, 1, sv_frame);
if stage_1_on_4_mics || stage_1_on_8_mics
    v3_frame = filter(b_air_v_3, 1, sv_frame);
    v4_frame = filter(b_air_v_4, 1, sv_frame);
    v5_frame = filter(b_air_v_5, 1, sv_frame);
    v6_frame = filter(b_air_v_6, 1, sv_frame);
    v7_frame = filter(b_air_v_7, 1, sv_frame);
    v8_frame = filter(b_air_v_8, 1, sv_frame);
end

% Add measurement noise
% NOTE: This is essentially regularizing (biasing) R_mc so dont make it unrealistically high
noise_mag = rms(sv_frame) .* 10 ^ (-1*SNR_dB / 20);
v1_frame = v1_frame + noise_mag .* randn(length(v1_frame), 1);
v2_frame = v2_frame + noise_mag .* randn(length(v2_frame), 1);
if stage_1_on_4_mics || stage_1_on_8_mics
    v3_frame = v3_frame + noise_mag .* randn(length(v3_frame), 1);
    v4_frame = v4_frame + noise_mag .* randn(length(v4_frame), 1);
    v5_frame = v5_frame + noise_mag .* randn(length(v5_frame), 1);
    v6_frame = v6_frame + noise_mag .* randn(length(v6_frame), 1);
    v7_frame = v7_frame + noise_mag .* randn(length(v7_frame), 1);
    v8_frame = v8_frame + noise_mag .* randn(length(v8_frame), 1);
end

% y1_frame = [y1_frame ; v1_frame];
% y2_frame = [y2_frame ; v2_frame];
% y3_frame = [y3_frame ; v3_frame];
% y4_frame = [y4_frame ; v4_frame];
% y5_frame = [y5_frame ; v5_frame];
% y6_frame = [y6_frame ; v6_frame];
% y7_frame = [y7_frame ; v7_frame];
% y8_frame = [y8_frame ; v8_frame];
% 
% s_frame = [s_frame ; sv_frame];


%% MINT: Ideal solution based on known AIRs

tau = 108;
h_air_list = [b_air_1 b_air_2];
MINT(h_air_list, tau);



%% LPC on clean speech (Not actually used in algo, just for reference)

L          = length(s_frame);
w_hamming  = hamming(L);
s_w = s_frame .* w_hamming;

% Zero padding (for autocorr method) -- not sure this is required but doesnt hurt
s_w = [zeros(p1,1) ; s_w ; zeros(p1,1)];

% Biased autocorr calc because unbiased xcorr introduces a changing scale factor (n-|m|)
% which differs from the set up of the autocorrelation method for LPC
% (required for a stable inverse filter)
[r_s, lags] = xcorr(s_w, s_w, 'biased', p1); 
idx_lag0    = find(lags==0);

R_s = toeplitz(r_s(idx_lag0:(idx_lag0+p1-1)));
r_s = r_s((idx_lag0+1):(idx_lag0+p1));

alpha_s = R_s \ r_s; % Option 1: Solve by gaussian elimination
%alpha_s = pinv(R_s) * r_s; % Option 2: Solve by pseudo-inverse

e_s = filter([1 ; -1*alpha_s], 1, s_frame);

[V_s,w_freqz] = freqz(1, [1 ; -1*alpha_s], Nfft); 
freqs_freqz   = w_freqz * (fs / (2*pi));


% LPC on reverberant speech (Solve Yule-walker via Gaussian elimination)
% y1_w        = y1_frame .* w_hamming;
% y1_w = [zeros(p1,1) ; y1_w ; zeros(p1,1)];
% 
% [phi_y1y1, lags] = xcorr(y1_w, y1_w, 'biased', p1);
% idx_lag0    = find(lags==0);
% 
% R_y1 = toeplitz(phi_y1y1(idx_lag0:(idx_lag0+p1-1)));
% r_y1 = phi_y1y1((idx_lag0+1):(idx_lag0+p1));
% 
% alpha_y1 = R_y1 \ r_y1;
% 
% e_y1 = filter([1 ; -1*alpha_y1], 1, y1_frame);

% 
% figure()
% subplot(3,1,1)
% plot(y1_frame)
% title('Frame of Speech')
% xlabel('samples')
% subplot(3,1,2)
% plot(e_y1)
% title('LPC Residual')
% xlabel('samples')
% subplot(3,1,3)
% plot(f_S_frame ./ 1000, 10*log10(Y1m_frame));
% hold on;
% plot(f_V_s ./ 1000,20*log10(abs(V_y1)));
% xlabel('Frequency [kHz]')
% ylabel('dB')
% title('Vocal Tract model')
% legend('Speech Spectrum', 'LPC Inverse Filter')
% sgtitle("E.g., Results of LPC on Reverberant Speech (One channel)")
 
figure()
subplot(3,1,1)
plot(s_frame)
title('Frame of Speech')
xlabel('samples')
subplot(3,1,2)
plot(e_s)
title('LPC Residual')
xlabel('samples')
subplot(3,1,3)
plot(freqs_welch ./ 1000, 10*log10(Sm_frame));
hold on;
plot(freqs_freqz ./ 1000, 20*log10(abs(V_s) .* (sqrt(max(Sm_frame)) / max(abs(V_s)))));
xlabel('Frequency [kHz]')
ylabel('dB')
title('Vocal Tract model')
legend('Clean Speech Spectrum', 'LPC Inverse Filter')
sgtitle("EXAMPLE ONLY: Results of LPC on Clean Speech (Source, not reverberant)")

% figure()
% zplane(1, [1 ; -1*alpha_s]')
% title('LPC inverse filter (clean speech)')
% 
% figure()
% zplane(1, [1 ; -1*alpha_y1]')
% title('LPC inverse filter (reverberant speech)')

%% Stage 1: Source Whitening

L          = length(y1_frame);
w_hamming  = hamming(L);

y1_w = y1_frame .* w_hamming;
y2_w = y2_frame .* w_hamming;
if stage_1_on_4_mics || stage_1_on_8_mics
    y3_w = y3_frame .* w_hamming;
    y4_w = y4_frame .* w_hamming;
    y5_w = y5_frame .* w_hamming;
    y6_w = y6_frame .* w_hamming;
    y7_w = y7_frame .* w_hamming;
    y8_w = y8_frame .* w_hamming;
end

y1_w = [zeros(p1,1) ; y1_w ; zeros(p1,1)];
y2_w = [zeros(p1,1) ; y2_w ; zeros(p1,1)];
if stage_1_on_4_mics || stage_1_on_8_mics
    y3_w = [zeros(p1,1) ; y3_w ; zeros(p1,1)];
    y4_w = [zeros(p1,1) ; y4_w ; zeros(p1,1)];
    y5_w = [zeros(p1,1) ; y5_w ; zeros(p1,1)];
    y6_w = [zeros(p1,1) ; y6_w ; zeros(p1,1)];
    y7_w = [zeros(p1,1) ; y7_w ; zeros(p1,1)];
    y8_w = [zeros(p1,1) ; y8_w ; zeros(p1,1)];
end

[phi_y1y1, lags] = xcorr(y1_w, y1_w, 'biased', p1);
[phi_y2y2,    ~] = xcorr(y2_w, y2_w, 'biased', p1);
if stage_1_on_4_mics || stage_1_on_8_mics
    [phi_y3y3,    ~] = xcorr(y3_w, y3_w, 'biased', p1);
    [phi_y4y4,    ~] = xcorr(y4_w, y4_w, 'biased', p1);
    [phi_y5y5,    ~] = xcorr(y5_w, y5_w, 'biased', p1);
    [phi_y6y6,    ~] = xcorr(y6_w, y6_w, 'biased', p1);
    [phi_y7y7,    ~] = xcorr(y7_w, y7_w, 'biased', p1);
    [phi_y8y8,    ~] = xcorr(y8_w, y8_w, 'biased', p1);
end
idx_lag0     = find(lags==0);

if stage_1_on_8_mics
    phi_ym = phi_y1y1 + phi_y2y2 + phi_y3y3 + phi_y4y4 + phi_y5y5 + phi_y6y6 + phi_y7y7 + phi_y8y8;
elseif stage_1_on_4_mics
    phi_ym = phi_y1y1 + phi_y2y2 + phi_y3y3 + phi_y4y4;
else
    phi_ym = phi_y1y1 + phi_y2y2;
end

R_ym = toeplitz(phi_ym(idx_lag0:(idx_lag0+p1-1)));
r_ym = phi_ym((idx_lag0+1):(idx_lag0+p1));

alpha_ym = R_ym \ r_ym; % Option 1: Solve by Gaussian elimination
%alpha_ym = pinv(R_ym) * r_ym; % Option 2: Solve by pseudo-inverse

e_y1_m = filter([1 ; -1*alpha_ym], 1, y1_frame);

% Compute source-whitened reverberant microphone signals
if (stage_1_on_clean_speech)
    x1_frame = filter([1 ; -1*alpha_s], 1, y1_frame);
    x2_frame = filter([1 ; -1*alpha_s], 1, y2_frame);
else
    x1_frame = filter([1 ; -1*alpha_ym], 1, y1_frame);
    x2_frame = filter([1 ; -1*alpha_ym], 1, y2_frame);
end

%% Stage 1 Plots

Y1m_frame           = pwelch(y1_frame, Nwin_welch, Nover_welch, Nfft_welch);
X1                  = pwelch(x1_frame, Nwin_welch, Nover_welch, Nfft_welch);
X2                  = pwelch(x2_frame, Nwin_welch, Nover_welch, Nfft_welch);

[V_ym,~] = freqz(1, [1 ; -1*alpha_ym], Nfft); 

figure()
subplot(3,1,1)
plot(y1_frame)
title('Frame of Speech')
xlabel('samples')
subplot(3,1,2)
plot(e_y1_m)
title('LPC Residual')
xlabel('samples')
subplot(3,1,3)
plot(freqs_welch ./ 1000, 10*log10(Sm_frame));
hold on;
plot(freqs_freqz ./ 1000,20*log10(abs(V_ym .* (sqrt(max(Sm_frame)) / max(abs(V_ym))))));
xlabel('Frequency [kHz]')
ylabel('dB')
title('Vocal Tract model')
legend('Clean Speech Spectrum (Not reverberant)', 'LPC Inverse Filter')
sgtitle("Stage 1: Results of LPC on Reverberant Speech (Avg over channels)")

figure()
subplot(2,1,1)
plot(freqs_welch, 10*log10(X1));
ylabel('dB')
xlabel('Frequency [Hz]')
title('Source Whitened Reverberant Signal 1 Spectrum (X1)')
subplot(2,1,2)
plot(freqs_welch, 10*log10(X2));
ylabel('dB')
xlabel('Frequency [Hz]')
title('Source Whitened Reverberant Signal 2 Spectrum (X2)')

%% STAGE 2: MULTI-CHANNEL LINEAR PREDICTION ON SOURCE-WHITENED REVERBERANT SIGNALS

if test_stage2_alone
    % Perform MC-LPC directly on reverberant signals (not source-whitened)
    x1_frame = y1_frame;
    x2_frame = y2_frame;
end

h0 = zeros(2,1); % First vector coefficient (estimate) of multi-path channel -- for time alignment
delay_comp = -1;
delay_filter_1 = [1];
delay_filter_2 = [1];
coeff_thresh = 0.1;

while min(abs(h0)) < coeff_thresh

    delay_comp = delay_comp + 1;
    if (abs(h0(1)) > coeff_thresh)
        x1_frame = [zeros(1,1) ; x1_frame(1:(end-1))];
        delay_filter_1 = [0; delay_filter_1];
    end
    if (abs(h0(2)) > coeff_thresh)
        x2_frame = [zeros(1,1) ; x2_frame(1:(end-1))];
        delay_filter_2 = [0; delay_filter_2];
    end    
    
    L          = length(x1_frame);
    w_hamming  = hamming(L);

    x1_w = x1_frame .* w_hamming;
    x2_w = x2_frame .* w_hamming;
    
    x1_w = [zeros(p2,1) ; x1_w ; zeros(p2,1)];
    x2_w = [zeros(p2,1) ; x2_w ; zeros(p2,1)];
    
    [phi_x1x1, lags] = xcorr(x1_w, x1_w, 'biased', p2);
    [phi_x1x2,    ~] = xcorr(x1_w, x2_w, 'biased', p2);
    
    [phi_x2x1,    ~] = xcorr(x2_w, x1_w, 'biased', p2);
    [phi_x2x2,    ~] = xcorr(x2_w, x2_w, 'biased', p2);
    
    idx_lag0       = find(lags==0);
    
    R_mc = zeros(2*p2, 2*p2);
    r_mc = zeros(2,    2*p2);
    
    k = 0;
    for block_row = 1:p2
        
        for block_col = 1:p2
            R_xx_k = [phi_x1x1(idx_lag0+k) phi_x1x2(idx_lag0+k) ; ...
                      phi_x2x1(idx_lag0+k) phi_x2x2(idx_lag0+k)];
    
            row_0 = (block_row - 1) * 2 + 1;
            row_1 = (block_row) * 2;
    
            col_0 = (block_col - 1) * 2 + 1;
            col_1 = (block_col) * 2;
    
            R_mc((row_0:row_1), (col_0:col_1)) = R_xx_k;
    
            k = k + 1;
        end
    
        k = k - p2 - 1 ;
    end
    
    k = 1;
    for block_col = 1:p2
        R_xx_k = [phi_x1x1(idx_lag0+k) phi_x1x2(idx_lag0+k) ; ...
                  phi_x2x1(idx_lag0+k) phi_x2x2(idx_lag0+k)];
    
        col_0 = (block_col - 1) * 2 + 1;
        col_1 = (block_col) * 2;
    
        r_mc(:, col_0:col_1) = R_xx_k;
    
        k = k+1;
    end
    
    % Solve Normal Equations
    % Note row formulation of system is used here, so aR = r (rather than Ra = r)
    %alpha_mc = r_mc / R_mc; % Option: Solve by Gaussian Elimination
    alpha_mc = r_mc * pinv(R_mc); % Option 2: Solve by pseudo inverse
    %alpha_mc = r_mc * ((inv(R_mc'*R_mc))*R_mc'); % Pseudo inverse for over-determined system
    %alpha_mc = r_mc * (R_mc'*(inv(R_mc*R_mc'))); % Pseudo inverse for under-determined system
    
    % Estimator Filters
    h_pred_1_from_1 = [0 alpha_mc(1, 1:2:(2*p2-1))];
    h_pred_1_from_2 = [0 alpha_mc(1, 2:2:2*p2)];
    
    h_pred_2_from_1 = [0 alpha_mc(2, 1:2:(2*p2-1))];
    h_pred_2_from_2 = [0 alpha_mc(2, 2:2:2*p2)];

    M = 2;
    h_pred_1_from_1 = [0 alpha_mc(1, 1:M:((p2-1)*M+1))];
    h_pred_1_from_2 = [0 alpha_mc(1, 2:M:((p2-1)*M+2))];
    % h_pred_1_from_3 = [0 alpha_mc(1, 3:M:((M-1)*p2+3))];
    % h_pred_1_from_4 = [0 alpha_mc(1, 4:M:((M-1)*p2+4))];
    % h_pred_1_from_5 = [0 alpha_mc(1, 5:M:((M-1)*p2+5))];
    % h_pred_1_from_6 = [0 alpha_mc(1, 6:M:((M-1)*p2+6))];

    h_pred_2_from_1 = [0 alpha_mc(2, 1:M:((p2-1)*M+1))];
    h_pred_2_from_2 = [0 alpha_mc(2, 2:M:((p2-1)*M+2))];

    x1_est = filter(h_pred_1_from_1, 1, x1_frame) + ...
             filter(h_pred_1_from_2, 1, x2_frame);
    
    x2_est = filter(h_pred_2_from_1, 1, x1_frame) + ...
             filter(h_pred_2_from_2, 1, x2_frame);
    
    e_x1_est = x1_frame - x1_est;
    e_x2_est = x2_frame - x2_est;
    
    % Estimate first vector coefficient of SIMO channel (for linear combiner)
    L          = length(e_x1_est);
    w_hamming  = hamming(L);
    x1_res_w   = e_x1_est .* w_hamming;
    x2_res_w   = e_x2_est .* w_hamming;
    
    [phi_x1x1_res, lags] = xcorr(x1_res_w, x1_res_w, 'biased', p2);
    [phi_x1x2_res,    ~] = xcorr(x1_res_w, x2_res_w, 'biased', p2);
    [phi_x2x1_res,    ~] = xcorr(x2_res_w, x1_res_w, 'biased', p2);
    [phi_x2x2_res,    ~] = xcorr(x2_res_w, x2_res_w, 'biased', p2);
    idx_lag0 = find(lags==0);
    
    R_xx_res_0 = [phi_x1x1_res(idx_lag0) phi_x1x2_res(idx_lag0) ; ...
                  phi_x2x1_res(idx_lag0) phi_x2x2_res(idx_lag0)];
    
    [V,D] = eig(R_xx_res_0);
    h0 = V(:, find(diag(D) == max(diag(D))));

    if apply_time_alignment ~= true
        break;
    end
end

if apply_time_alignment

    diff = abs(length(delay_filter_1) - length(delay_filter_2));
    if length(delay_filter_1) > length(delay_filter_2)
        delay_filter_2 = [delay_filter_2 ; zeros(diff,1)];
    else
        delay_filter_1 = [delay_filter_1 ; zeros(diff,1)];
    end

    % Adjust other signals for delay compensation
    s_frame  = filter([zeros(delay_comp,1) ; 1], 1, s_frame);
    y1_frame = filter(delay_filter_1, 1, y1_frame);
    y2_frame = filter(delay_filter_2, 1, y2_frame);
end


%% Stage 3: Inverse Filtering

% Compute DAP output (de-reverberate)
if test_stage2_alone
    s_est_1_dap = e_x1_est;
    s_est_2_dap = e_x2_est; % Only equal for white source (otherwise need to re-apply filter to non-whitened reverberant signal)
    s_est_dap = s_est_1_dap .* h0(1) + s_est_2_dap .* h0(2);

else
    s_est_1_dap = y1_frame - ...
              ((filter(h_pred_1_from_1, 1, y1_frame) + ...
              filter(h_pred_1_from_2, 1, y2_frame)));
    
    s_est_2_dap = y2_frame - ...
              ((filter(h_pred_2_from_1, 1, y1_frame) + ...
              filter(h_pred_2_from_2, 1, y2_frame)));
    
    s_est_dap = s_est_1_dap .* h0(1) + s_est_2_dap .* h0(2);

end

% Manually compensate gain ambiguity problem for error calc
if adjust_gain_ambiguity
    s_est_dap = s_est_dap .* (rms(s_frame) / rms(s_est_dap));
end

e_est_dap = s_frame - s_est_dap;

Sm_frame          = pwelch(s_frame,   Nwin_welch, Nover_welch, Nfft_welch);
S_est_dap_m_frame = pwelch(s_est_dap, Nwin_welch, Nover_welch, Nfft_welch);
E_est_dap_m_frame = pwelch(e_est_dap, Nwin_welch, Nover_welch, Nfft_welch);

ylim_max = max([max(abs(s_frame)), max(abs(s_est_dap)), max(abs(e_est_dap))]) * 1.2;
figure()
subplot(3,1,1)
plot(s_frame)
ylim([-1*ylim_max, ylim_max])
title('Source Signal')
subplot(3,1,2)
plot(s_est_dap)
ylim([-1*ylim_max, ylim_max])
title('Estimated Source Signal (Delay-and-Predict)')
subplot(3,1,3)
plot(e_est_dap)
%ylim([-1*ylim_max, ylim_max])
title('Error')
sgtitle('Stage 2: MC-LPC Results')

figure()
subplot(2,1,1)
plot(freqs_welch ./ 1000, 10*log10(Sm_frame))
hold on;
plot(freqs_welch ./ 1000, 10*log10(S_est_dap_m_frame))
legend('Source', 'Estimated Source')
title('Magnitude Spectrum')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(freqs_welch ./ 1000, 10*log10(E_est_dap_m_frame))
title('Error Spectrum')
xlabel('Frequency [kHz]')
ylabel('dB')
sgtitle('Stage 2: MC-LPC Results')


%% Generate/Plot IR for total DAP system (pass impulse through)

s_impulse = zeros(L_channel + p1 + p2,1);
s_impulse(1) = 1;

h_channel_1 = filter(b_air_1, a_air_1, s_impulse);
h_channel_2 = filter(b_air_2, a_air_2, s_impulse);

x1_frame = filter(b_air_1, a_air_1, s_impulse);
x2_frame = filter(b_air_2, a_air_2, s_impulse);

x1_frame = filter(delay_filter_1, 1, x1_frame);
x2_frame = filter(delay_filter_2, 1, x2_frame);

h_1_dap = x1_frame - ...
          ((filter(h_pred_1_from_1, 1, x1_frame) + ...
          filter(h_pred_1_from_2, 1, x2_frame)));

h_2_dap = x2_frame - ...
          ((filter(h_pred_2_from_1, 1, x1_frame) + ...
          filter(h_pred_2_from_2, 1, x2_frame)));

h_dap = h_1_dap .* h0(1) + h_2_dap .* h0(2);

H_dap   = fft(h_dap, Nfft);
Hm_dap  = abs(H_dap(1:(Nfft/2)));
Hp_dap  = angle(H_dap(1:(Nfft/2)));
freqs_fft = (0:(length(Hm_dap)-1))' .* (fs / Nfft);

H_channel_1  = freqz(b_air_1, a_air_1, Nfft);
Hm_channel_1 = abs(H_channel_1);
Hp_channel_1  = angle(H_channel_1);

H_channel_2  = freqz(b_air_2, a_air_2, Nfft);
Hm_channel_2 = abs(H_channel_2);
Hp_channel_2  = angle(H_channel_2);

H_pred_1_from_1  = freqz(h_pred_1_from_1, 1, Nfft);
Hm_pred_1_from_1 = abs(H_pred_1_from_1);
Hp_pred_1_from_1 = angle(H_pred_1_from_1);

H_pred_1_from_2  = freqz(h_pred_1_from_2, 1, Nfft);
Hm_pred_1_from_2 = abs(H_pred_1_from_2);
Hp_pred_1_from_2 = angle(H_pred_1_from_2);

H_pred_2_from_1  = freqz(h_pred_2_from_1, 1, Nfft);
Hm_pred_2_from_1 = abs(H_pred_2_from_1);
Hp_pred_2_from_1 = angle(H_pred_2_from_1);

H_pred_2_from_2  = freqz(h_pred_2_from_2, 1, Nfft);
Hm_pred_2_from_2 = abs(H_pred_2_from_2);
Hp_pred_2_from_2 = angle(H_pred_2_from_2);

figure()
subplot(5,1,1)
plot((0:(length(h_channel_1)-1)) .* (1/fs), h_channel_1);
xlabel('Time [sec]')
title("Channel 1 Impulse Response")
subplot(5,1,2)
plot((0:(length(h_channel_2)-1)) .* (1/fs), h_channel_2);
xlabel('Time [sec]')
title("Channel 2 Impulse Response")
subplot(5,1,3)
plot((0:(length(h_1_dap)-1)) .* (1/fs), h_1_dap);
xlabel('Time [sec]')
title("h 1 dap Equalized Impulse Response")
subplot(5,1,4)
plot((0:(length(h_2_dap)-1)) .* (1/fs), h_2_dap);
xlabel('Time [sec]')
title("h 2 dap Equalized Impulse Response")
subplot(5,1,5)
plot((0:(length(h_dap)-1)) .* (1/fs), h_dap);
xlabel('Time [sec]')
title("h dap Equalized Impulse Response")

figure()
subplot(2,1,1)
stem((0:(length(h_channel_1)-1)) .* (1/fs), h_channel_1);
xlabel('Time [sec]')
title("Channel 1 Impulse Response")
subplot(2,1,2)
stem((0:(length(h_dap)-1)) .* (1/fs), h_dap);
xlabel('Time [sec]')
title("Equalized Impulse Response (Stem Plot)")

figure()
subplot(2,1,1)
plot((0:(length(h_channel_1)-1)) .* (1/fs), h_channel_1);
xlabel('Time [sec]')
title("Channel 1 Impulse Response")
subplot(2,1,2)
plot((0:(length(h_dap)-1)) .* (1/fs), h_dap);
xlabel('Time [sec]')
title("Equalized Impulse Response")


figure()
subplot(2,2,1)
plot(freqs_freqz ./ 1000, 20*log10(Hm_channel_1));
title('Channel 1 magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,3)
plot(freqs_freqz ./ 1000, Hp_channel_1);
title('Channel 1 phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

subplot(2,2,2)
plot(freqs_freqz ./ 1000, 20*log10(Hm_channel_2));
title('Channel 2 magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,4)
plot(freqs_freqz ./ 1000, Hp_channel_1);
title('Channel 2 phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

figure()
subplot(2,2,1)
plot(freqs_freqz ./ 1000, 20*log10(Hm_pred_1_from_1))
title('Prediction Filter h pred 1 from 1')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,3)
plot(freqs_freqz ./ 1000, 20*log10(Hm_pred_1_from_2))
title('Prediction Filter h pred 1 from 2')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,2)
plot(freqs_freqz ./ 1000, 20*log10(Hm_pred_2_from_1))
title('Prediction Filter h pred 2 from 1')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,4)
plot(freqs_freqz ./ 1000, 20*log10(Hm_pred_2_from_2))
title('Prediction Filter h pred 2 from 2')
xlabel('Frequency [kHz]')
ylabel('dB')


% figure()
% subplot(1,2,1)
% zplane(b_air_1', a_air_1')
% title('Channel 1')
% subplot(1,2,2)
% zplane(b_air_2', a_air_2')
% title('Channel 2')

figure()
subplot(2,1,1)
plot(freqs_fft ./ 1000, 20*log10(Hm_dap));
title('Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(freqs_fft ./ 1000, unwrap(Hp_dap));
title('Equalized phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

%% Energy Decay Curve Comparison

energy_air_1 = h_channel_1' * h_channel_1;
energy_dap   = h_dap' * h_dap;
edc_air_1    = zeros(length(h_channel_1), 1);
edc_dap      = zeros(length(h_channel_1), 1);
for n = 1:length(h_channel_1)

    edc_air_1(n) = energy_air_1;
    edc_dap(n)   = energy_dap;

    energy_air_1 = energy_air_1 - h_channel_1(n)^2;
    energy_dap   = energy_dap - h_dap(n)^2;
end

figure()
plot((0:(length(h_channel_1)-1)) .* (1/fs), 20*log10(edc_air_1));
hold on;
plot((0:(length(h_channel_1)-1)) .* (1/fs), 20*log10(edc_dap .* (max(edc_air_1) / max(edc_dap))));
ylim([-80 60])
%xlim([0 (length(h_channel_1) * (1/fs))])
xlabel('Time [sec]')
ylabel('Energy remaining [dB]')
legend('Reverb Energy', 'Equalized Reverb Energy')
title('Energy Decay Curve')

%% Temp

[s_frame,fs] = audioread("SA1.wav");
y1_frame = filter(b_air_1, a_air_1, s_frame);
y2_frame = filter(b_air_2, a_air_2, s_frame);
s_est_1_dap = y1_frame - ...
          ((filter(h_pred_1_from_1, 1, y1_frame) + ...
          filter(h_pred_1_from_2, 1, y2_frame)));

s_est_2_dap = y2_frame - ...
          ((filter(h_pred_2_from_1, 1, y1_frame) + ...
          filter(h_pred_2_from_2, 1, y2_frame)));

s_est_dap = s_est_1_dap .* h0(1) + s_est_2_dap .* h0(2);

