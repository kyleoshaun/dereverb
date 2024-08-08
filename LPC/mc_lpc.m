addpath('../../samples/')
addpath('../../AIR_Databases/AIR_1_4_BinauralDatabase/')

close all
clear
clc

%rng("shuffle");
%rng("default");
rng(7)

%% Parameters

fs = 16000;
p1 = 100; % Source Whitening Order
p2 = 20; % MC-LPC Order
test_stage2_alone = false;
adjust_gain_ambiguity = true;
SNR_dB = 120;

Nfft = 4096;


[s_frame,fs] = audioread("assgn1.wav");
%[s_frame,fs] = audioread("SA1_phoneme.wav");
%[s_frame,fs] = audioread("SA1.wav");
%s_frame = s_frame((1.6*fs):(3.05*fs));
if test_stage2_alone
    % Set s_frame to white noise (rather than speech) to debug MC-LPC alone
    s_frame = randn(500000,1);
end

% exponential decay curve
%L_channel = 10;
%tau = L_channel/4;
%exp_decay = exp(-1 .* (1:L_channel)' ./ tau);

% Channel
b_air_1 = [1 -0.50  0.42  0.05 -0.44     0]';
b_air_2 = [1 -0.10  0.40 -0.40  0.10  0.05]';
a_air_1 = 1; %[1 -0.50  0.42  0.05 -0.44     0]';
a_air_2 = 1; %[1 -0.10  0.40 -0.40  0.10  0.05]';
L_channel = length(b_air_1);

% Real Binaural AIR
% [b_air_1, b_air_2] = load_real_AIR(fs);
% a_air_1 = 1;
% a_air_2 = 1;

%Exponential shortening of AIR (To simplify)
%b_air_1 = b_air_1(1:L_channel) .* exp_decay;
%b_air_2 = b_air_2(1:L_channel) .* exp_decay;


%b_air_1 = [flip(b_air_1) ; b_air_1];
%b_air_2 = [flip(b_air_2) ; b_air_2];
%b_air_1 = [zeros(2000,1) ; b_air_1];
%b_air_2 = [zeros(2000,1) ; b_air_2];
L_channel = length(b_air_1);

% Compute reverberant signals
y1_frame = filter(b_air_1, a_air_1, s_frame);
y2_frame = filter(b_air_2, a_air_2, s_frame);

% Add measurement noise
% NOTE: This is essentially regularizing (biasing) R_mc so dont make it unrealistically high
noise_mag = rms(s_frame) .* 10 ^ (-1*SNR_dB / 20);
y1_frame = y1_frame + noise_mag .* randn(length(y1_frame), 1);
y2_frame = y2_frame + noise_mag .* randn(length(y2_frame), 1);

%% STAGE 1: SOURCE WHITENING

L          = length(s_frame);
w_hamming  = hamming(L);

% LPC on clean speech (Solve Yule-walker via Gaussian elimination)
s_w = s_frame .* w_hamming;

% Zero padding (for autocorr method) -- not sure this is required but doesnt hurt
s_w = [zeros(p1,1) ; s_w ; zeros(p1,1)];

% Biased autocorr calc because unbiased xcorr introduces a changing scale factor (n-|m|)
% which differs from the set up of the autocorrelation method for LPC
% (required for a stable inverse filter)
[r_s, lags] = xcorr(s_w, s_w, 'biased'); 
idx_lag0    = find(lags==0);

R_s = toeplitz(r_s(idx_lag0:(idx_lag0+p1-1)));
r_s = r_s((idx_lag0+1):(idx_lag0+p1));

alpha_s = R_s \ r_s;

e_s = filter([1 ; -1*alpha_s], 1, s_frame);

% LPC on reverberant speech (Solve Yule-walker via Gaussian elimination)
y1_w        = y1_frame .* w_hamming;
y1_w = [zeros(p1,1) ; y1_w ; zeros(p1,1)];

[phi_y1y1, lags] = xcorr(y1_w, y1_w, 'biased');
idx_lag0    = find(lags==0);

R_y1 = toeplitz(phi_y1y1(idx_lag0:(idx_lag0+p1-1)));
r_y1 = phi_y1y1((idx_lag0+1):(idx_lag0+p1));

alpha_y1 = R_y1 \ r_y1;

e_y1 = filter([1 ; -1*alpha_y1], 1, y1_frame);

% Multi-channel LPC on reverberant speech (Solve Yule-walker via Gaussian elimination)
y1_w        = y1_frame .* w_hamming;
y2_w        = y2_frame .* w_hamming;

y1_w = [zeros(p1,1) ; y1_w ; zeros(p1,1)];
y2_w = [zeros(p1,1) ; y2_w ; zeros(p1,1)];

[phi_y1y1, lags] = xcorr(y1_w, y1_w, 'biased');
[phi_y2y2,    ~] = xcorr(y2_w, y2_w, 'biased');
idx_lag0     = find(lags==0);

phi_y_avg = phi_y1y1 + phi_y2y2;

R_y_avg = toeplitz(phi_y_avg(idx_lag0:(idx_lag0+p1-1)));
r_y_avg = phi_y_avg((idx_lag0+1):(idx_lag0+p1));

alpha_y_avg = R_y_avg \ r_y_avg;

e_y1_avg = filter([1 ; -1*alpha_y_avg], 1, y1_frame);


% Compute source-whitened reverberant microphone signals
x1_frame = filter([1 ; -1*alpha_y_avg], 1, y1_frame);
x2_frame = filter([1 ; -1*alpha_y_avg], 1, y2_frame);

% Plots
[V_s,w_freqz] = freqz(1, [1 ; -1*alpha_s], Nfft); 
freqs_freqz   = w_freqz * (fs / (2*pi));

[X1, ~] = pwelch(x1_frame, 2048, 512, Nfft);
[X2, ~] = pwelch(x2_frame, 2048, 512, Nfft);
[Sm_frame, w_welch] = pwelch(s_frame, 2048, 512, Nfft);
Y1m_frame           = pwelch(y1_frame, 2048, 512, Nfft);
freqs_welch         = w_welch .* (fs / (2*pi));


[V_y1, ~] = freqz(1, [1 ; -1*alpha_y1], Nfft); 
[V_ym, ~] = freqz(1, [1 ; -1*alpha_y_avg], Nfft); 

plot_scale = max(Sm_frame) / max(abs(V_s));
 
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
plot(freqs_welch, 10*log10(Sm_frame));
hold on;
plot(freqs_freqz,20*log10(abs(V_s)));
xlabel('Frequency [Hz]')
ylabel('dB')
title('Vocal Tract model')
legend('Clean Speech Spectrum', 'LPC Inverse Filter')
sgtitle("E.g., Results of LPC on Clean Speech (Source, not reverberant)")

figure()
subplot(3,1,1)
plot(y1_frame)
title('Frame of Speech')
xlabel('samples')
subplot(3,1,2)
plot(e_y1)
title('LPC Residual')
xlabel('samples')
subplot(3,1,3)
plot(freqs_welch, 10*log10(Y1m_frame));
hold on;
plot(freqs_freqz,20*log10(abs(V_y1)));
xlabel('Frequency [Hz]')
ylabel('dB')
title('Vocal Tract model')
legend('Speech Spectrum', 'LPC Inverse Filter')
sgtitle("E.g., Results of LPC on Reverberant Speech (One channel)")

figure()
subplot(3,1,1)
plot(y1_frame)
title('Frame of Speech')
xlabel('samples')
subplot(3,1,2)
plot(e_y1_avg)
title('LPC Residual')
xlabel('samples')
subplot(3,1,3)
plot(freqs_welch, 10*log10(Y1m_frame));
hold on;
plot(freqs_freqz,20*log10(abs(V_ym)));
xlabel('Frequency [Hz]')
ylabel('dB')
title('Vocal Tract model')
legend('Speech Spectrum (channel 1)', 'LPC Inverse Filter')
sgtitle("Stage 1: Results of LPC on Reverberant Speech (Avg over channels)")

figure()
subplot(2,1,1)
plot(freqs_welch, 20*log10(X1));
xlabel('Frequency [Hz]')
xlabel('dB')
title('Source Whitened Reverberant Signal 1 (x1)')
subplot(2,1,2)
plot(freqs_welch, 20*log10(X2));
xlabel('Frequency [Hz]')
xlabel('dB')
title('Source Whitened Reverberant Signal 1 (x2)')


% figure()
% zplane(1, [1 ; -1*alpha_s]')
% title('LPC inverse filter (clean speech)')
% 
% figure()
% zplane(1, [1 ; -1*alpha_y1]')
% title('LPC inverse filter (reverberant speech)')
% 
% figure()
% zplane(1, [1 ; -1*alpha_ym]')
% title('Multi Channel LPC inverse filter (reverberant speech)')

%% STAGE 2: MULTI-CHANNEL LINEAR PREDICTION ON SOURCE-WHITENED REVERBERANT SIGNALS

if test_stage2_alone
    % Debug Tool: perform MC-LPC directly on reverberant signals (not source-whitened)
    x1_frame = y1_frame;
    x2_frame = y2_frame;
end

L          = length(s_frame);
w_hamming  = hamming(L);
x1_w       = x1_frame .* w_hamming;
x2_w       = x2_frame .* w_hamming;

x1_w = [zeros(p2,1) ; x1_w ; zeros(p2,1)];
x2_w = [zeros(p2,1) ; x2_w ; zeros(p2,1)];

% x_mc = E[x1(n-1) ... x1(n-p2) x2(n) ... x2(n-p2) ... xM(n-p2)]'
% r_mc = E[x_mc * d(n)] = E[x_mc * x1(n)]
%      = E[x1(n-1)x1(n) ... x1(n-p2)x1(n) x2(n)x1(n) ... x2(n-p2)d(n) ... xM(n-p2)d(n)]'
% R_mc = E[x_mc * x_mc']
%         [x1(n-1)x1(n-1)  ... x1(n-1)x1(n-p2)  x1(n-1)x2(n)  ... x1(n-1)x2(n-p2)  ... x1(n-1)xM(n-p2)  ]
%      = E[x1(n-2)x1(n-1)  ... x1(n-2)x1(n-p2)  x1(n-2)x2(n)  ... x1(n-2)x2(n-p2)  ... x1(n-2)xM(n-p2)  ]
%         [      :                    :               :                 :                    :          ]
%         [xM(n-p2)x1(n-1) ... xM(n-p2)x1(n-p2) xM(n-p2)x2(n) ... xM(n-p2)x2(n-p2) ... xM(n-p2)xM(n-p2) ]

[phi_x1x1, lags] = xcorr(x1_w, x1_w, 'biased', p2);
[phi_x1x2,    ~] = xcorr(x1_w, x2_w, 'biased', p2);

[phi_x2x1,    ~] = xcorr(x2_w, x1_w, 'biased', p2);
[phi_x2x2,    ~] = xcorr(x2_w, x2_w, 'biased', p2);

idx_lag0       = find(lags==0);


R_x1x1 = toeplitz(phi_x1x1(idx_lag0:-1:(idx_lag0-p2+1)), phi_x1x1(idx_lag0:(idx_lag0+p2-1)));
R_x1x2 = toeplitz(phi_x1x2(idx_lag0:-1:(idx_lag0-p2+1)), phi_x1x2(idx_lag0:(idx_lag0+p2-1)));

R_x2x1 = toeplitz(phi_x2x1(idx_lag0:-1:(idx_lag0-p2+1)), phi_x2x1(idx_lag0:(idx_lag0+p2-1)));
R_x2x2 = toeplitz(phi_x2x2(idx_lag0:-1:(idx_lag0-p2+1)), phi_x2x2(idx_lag0:(idx_lag0+p2-1)));

% Enforcing Symmetry (MC-LP Thesis I found says this should be the case -- I checked and it was pretty close, error = 1e-16)
% R_x1x1 = toeplitz(phi_x1x1(idx_lag0:(idx_lag0+p2-1)));
% R_x2x2 = toeplitz(phi_x2x2(idx_lag0:(idx_lag0+p2-1)));
% 
% R_x1x2 = toeplitz(phi_x1x2(idx_lag0:-1:(idx_lag0-p2+1)), phi_x1x2(idx_lag0:(idx_lag0+p2-1)));
% R_x2x1 = R_x1x2';


% 2 Channel LPC
R_mc = [R_x1x1 R_x1x2 ; ...
        R_x2x1 R_x2x2];

r_mc = [phi_x1x1((idx_lag0-1):-1:(idx_lag0-p2))   ; ...
        phi_x2x1((idx_lag0-1):-1:(idx_lag0-p2))];

% Regularize
%psi = 0; %0.0001 * min(min(abs(R_mc)));
%fprintf("Condition number of Spatio-Temporal correlaton matrix = %.2e\n", cond(R_mc));
% R_mc = R_mc + psi*eye(size(R_mc));
% fprintf("Post-Regularization cond(R_mc + %.2e*I) = %.2e\n", psi, cond(R_mc));

alpha_mc = R_mc \ r_mc;
%alpha_mc = pinv(R_mc) * r_mc;

% Estimator Filters
eq_1 = [0 ; alpha_mc(((1-1)*p2+1):1*p2)];
eq_2 = [0 ; alpha_mc(((2-1)*p2+1):2*p2)];


x1_est = filter(eq_1, 1, x1_frame) + ...
         filter(eq_2, 1, x2_frame);

e_x1_est = x1_frame - x1_est;

%% Stage 3: Inverse Filtering

% Compute DAP output (de-reverberate)
if test_stage2_alone
    s_est_dap = e_x1_est; % Only equal for white source (otherwise need to re-apply filter to non-whitened reverberant signal)
else
s_est_dap  = y1_frame - ...
             (filter(eq_1, 1, y1_frame) + ...
              filter(eq_2, 1, y2_frame));
end



% Manually compensate gain ambiguity problem for error calc
if adjust_gain_ambiguity
    s_est_dap = s_est_dap .* (max(s_frame) / max(s_est_dap));
end

%e_est_dap = s_frame - s_est_dap;

Sm_frame          = pwelch(s_frame,   512, 256, Nfft);
S_est_dap_m_frame = pwelch(s_est_dap, 512, 256, Nfft);
%E_est_dap_m_frame = pwelch(e_est_dap, 512, 256, Nfft);

figure()
subplot(2,1,1)
plot((0:(Nfft/2)) .* (fs/Nfft), 10*log10(Sm_frame))
hold on;
plot((0:(Nfft/2)) .* (fs/Nfft), 10*log10(S_est_dap_m_frame))
xlabel('Frequency [Hz]')
legend('Input (Gain-Compensated', 'Equalized Input')
subplot(2,1,2)
plot((0:(Nfft/2)) .* (fs/Nfft),  10*log10(S_est_dap_m_frame) - 10*log10(Sm_frame))
legend('Error (PSD EQd input - PSD input)')
%ylim([-1, 1])
sgtitle('Equalization Results')

%% Generate/Plot IR for total DAP system (pass impulse through)

%s_impulse = zeros(round(L_channel*5),1);
s_impulse = zeros(round(L_channel+p2),1);
s_impulse(1) = 1;

h_channel_1 = filter(b_air_1, a_air_1, s_impulse);
h_channel_2 = filter(b_air_2, a_air_2, s_impulse);

x1_frame = filter(b_air_1, a_air_1, s_impulse);
x2_frame = filter(b_air_2, a_air_2, s_impulse);
x1_est = filter(eq_1, 1, x1_frame) + ...
         filter(eq_2, 1, x2_frame);

h_dap = x1_frame - x1_est;

H_dap   = fft(h_dap, Nfft);
Hm_dap  = abs(H_dap(1:(Nfft/2)));
Hp_dap  = angle(H_dap(1:(Nfft/2)));
f_H_air = (0:(length(Hm_dap)-1))' .* (fs / Nfft);

H_channel_1  = freqz(b_air_1, a_air_1, Nfft, 'whole');
Hm_channel_1 = abs(H_channel_1(1:(Nfft/2)));
Hp_channel_1  = angle(H_channel_1(1:(Nfft/2)));

H_channel_2  = freqz(b_air_2, a_air_2, Nfft, 'whole');
Hm_channel_2 = abs(H_channel_2(1:(Nfft/2)));
Hp_channel_2  = angle(H_channel_2(1:(Nfft/2)));

% 
% figure()
% subplot(2,1,1)
% stem((0:(length(h_channel_1)-1)) .* (1/fs), h_channel_1);
% xlabel('Time [sec]')
% title("Channel 1 Impulse Response")
% subplot(2,1,2)
% stem((0:(length(h_dap)-1)) .* (1/fs), h_dap);
% xlabel('Time [sec]')
% title("Equalized Impulse Response")

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
plot(f_H_air ./ 1000, 20*log10(Hm_channel_1));
title('Channel 1 magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,3)
plot(f_H_air ./ 1000, Hp_channel_1);
title('Channel 1 phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

subplot(2,2,2)
plot(f_H_air ./ 1000, 20*log10(Hm_channel_2));
title('Channel 2 magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,2,4)
plot(f_H_air ./ 1000, Hp_channel_1);
title('Channel 2 phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

% figure()
% subplot(1,2,1)
% zplane(b_air_1', a_air_1')
% title('Channel 1')
% subplot(1,2,2)
% zplane(b_air_2', a_air_2')
% title('Channel 2')

figure()
subplot(2,1,1)
plot(f_H_air ./ 1000, 20*log10(Hm_dap));
title('Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(f_H_air ./ 1000, unwrap(Hp_dap));
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


