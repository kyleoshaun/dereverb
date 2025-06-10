% This script completes multi-channel linear prediction on 
% two simulated reverberant microphone signals: x1 and x2
%
%     s(n) --> [ AIR 1 ] --> o| x1 
%  |< 
%     s(n) --> [ AIR 2 ] --> o| x2 
%
% Specifically, here I am performing two separate multi-channel linear predictions:
%   1. Predicting x1 from past samples of both x1 and x2
%   2. Predicting x2 from past samples of both x1 and x2
% Both of these linear prediction problems get packed into one matrix formulation
% and the solution to both are computed together. This is the same formulation
% as references I've found online E.g., https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6695272
%
% I also include the combination of the prediction error residuals from (1) and (2)
% using the "precursor coefficient" just like in the original Delay-and-Predict 
% dereverb paper from Dirk Slock https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1661221
%
% I formulate the normal equations by packing the input data as time-first,
% sensor # second. This is the opposite of the derivations I've found online.
% See mc_lpc_alt_packing.m for the opposite formulation.
%
% The multichannel linear prediction error filter (i.e., whitening filter)
% should remove the effect of the corresponding channel thus making it an equalizer
%
% To analyze the equalized response I generate an impulse and 
% pass it through the channel and then the prediction error filter
% 

addpath('../../samples/')
addpath('../../AIR_Databases/AIR_1_4_BinauralDatabase/')

close all
clear
clc

rng('shuffle')
seed = abs(round(randn(1,1) * 1000));
rng(724);

%% Parameters

p1 = 100; % Source Whitening Order
p2 = 200; % MC-LPC Order
test_stage2_alone = true;
stage_1_on_clean_speech = true;
adjust_gain_ambiguity = true;
SNR_dB = 300;


[s_frame,fs] = audioread("assgn1.wav");
s_frame = s_frame + 10^(-20/20) .* randn(length(s_frame),1);
%[s_frame,fs] = audioread("SA1_phoneme.wav");
%[s_frame,fs] = audioread("SA1.wav");
%s_frame = randn(500000,1);
if test_stage2_alone
    % Set s_frame to white noise (rather than speech) to debug MC-LPC alone
    s_frame = randn(500000,1);
end

Nfft = 4096;
Nfft_welch  = 2^(nextpow2(length(s_frame)));
Nfft_welch  = min(Nfft_welch, Nfft);
Nwin_welch  = Nfft_welch / 2;
Nover_welch = Nfft_welch / 4;

% exponential decay curve
L_channel = 20;
tau = L_channel/4;
exp_decay = exp(-1 .* (0:(L_channel-1))' ./ tau);

% Channel
b_air_2 = randn(L_channel,1);%[1 -1];%-0.50  0.42  0.05 -0.44     0]';
b_air_1 = randn(L_channel,1);%[1 1];%-0.10  0.40 -0.40  0.10  0.05]';
% b_air_3 = randn(L_channel,1);
% b_air_4 = randn(L_channel,1);
% b_air_5 = randn(L_channel,1);
% b_air_6 = randn(L_channel,1);
% b_air_7 = randn(L_channel,1);
% b_air_8 = randn(L_channel,1);
a_air_1 = 1; %[1 -0.50  0.42  0.05 -0.44     0]';
a_air_2 = 1; %[1 -0.10  0.40 -0.40  0.10  0.05]';
% a_air_3 = 1;
% a_air_4 = 1;
% a_air_5 = 1;
% a_air_6 = 1;
% a_air_7 = 1;
% a_air_8 = 1;

% Real Binaural AIR
% [b_air_1, b_air_2] = load_real_AIR(fs);
% a_air_1 = 1;
% a_air_2 = 1;

% Exponential shortening of AIR (To simplify)
b_air_1 = b_air_1(1:L_channel) .* exp_decay;
b_air_2 = b_air_2(1:L_channel) .* exp_decay;
% b_air_3 = b_air_3(1:L_channel) .* exp_decay;
% b_air_4 = b_air_4(1:L_channel) .* exp_decay;
% b_air_5 = b_air_5(1:L_channel) .* exp_decay;
% b_air_6 = b_air_6(1:L_channel) .* exp_decay;
% b_air_7 = b_air_7(1:L_channel) .* exp_decay;
% b_air_8 = b_air_8(1:L_channel) .* exp_decay;


%b_air_1 = [flip(b_air_1) ; b_air_1];
%b_air_2 = [flip(b_air_2) ; b_air_2];
%b_air_1 = [zeros(2000,1) ; b_air_1];
%b_air_2 = [zeros(2000,1) ; b_air_2];

% Add time delay
delay_1 = 0;
delay_2 = 100;
b_air_1 = [zeros(delay_1,1) ; b_air_1 ; zeros(delay_2,1)];
b_air_2 = [zeros(delay_2,1) ; b_air_2  ; zeros(delay_1,1)];

L_channel = length(b_air_1);

% Compute reverberant signals
x1_frame = filter(b_air_1, a_air_1, s_frame);
x2_frame = filter(b_air_2, a_air_2, s_frame);
% x3_frame = filter(b_air_3, a_air_3, s_frame);
% x4_frame = filter(b_air_4, a_air_4, s_frame);
% x5_frame = filter(b_air_5, a_air_5, s_frame);
% x6_frame = filter(b_air_6, a_air_6, s_frame);
% x7_frame = filter(b_air_7, a_air_7, s_frame);
% x8_frame = filter(b_air_8, a_air_8, s_frame);

% Add measurement noise
% NOTE: This is essentially regularizing (biasing) R_mc so dont make it unrealistically high
noise_mag = rms(s_frame) .* 10 ^ (-1*SNR_dB / 20);
x1_frame = x1_frame + noise_mag .* randn(length(x1_frame), 1);
x2_frame = x2_frame + noise_mag .* randn(length(x2_frame), 1);
% x3_frame = x3_frame + noise_mag .* randn(length(x3_frame), 1);
% x4_frame = x4_frame + noise_mag .* randn(length(x4_frame), 1);
% x5_frame = x5_frame + noise_mag .* randn(length(x5_frame), 1);
% x6_frame = x6_frame + noise_mag .* randn(length(x6_frame), 1);
% x7_frame = x7_frame + noise_mag .* randn(length(x7_frame), 1);
% x8_frame = x8_frame + noise_mag .* randn(length(x8_frame), 1);

%% Solution

L          = length(s_frame);
w_hamming  = hamming(L);
x1_w       = x1_frame .* w_hamming;
x2_w       = x2_frame .* w_hamming;

x1_w = [zeros(p2,1) ; x1_w ; zeros(p2,1)];
x2_w = [zeros(p2,1) ; x2_w ; zeros(p2,1)];

[phi_x1x1, lags] = xcorr(x1_w, x1_w, 'biased', p2);
[phi_x1x2,    ~] = xcorr(x1_w, x2_w, 'biased', p2);

[phi_x2x1,    ~] = xcorr(x2_w, x1_w, 'biased', p2);
[phi_x2x2,    ~] = xcorr(x2_w, x2_w, 'biased', p2);

idx_lag0 = find(lags==0);

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

alpha_mc = r_mc / R_mc; % Note row formulation of system is used here, so aR = r (rather than Ra = r)

% Estimator Filters
h_pred_1_from_1 = [0 alpha_mc(1, 1:2:(2*p2-1))];
h_pred_1_from_2 = [0 alpha_mc(1, 2:2:2*p2)];

h_pred_2_from_1 = [0 alpha_mc(2, 1:2:(2*p2-1))];
h_pred_2_from_2 = [0 alpha_mc(2, 2:2:2*p2)];



x1_est = filter(h_pred_1_from_1, 1, x1_frame) + ...
         filter(h_pred_1_from_2, 1, x2_frame);

x2_est = filter(h_pred_2_from_1, 1, x1_frame) + ...
         filter(h_pred_2_from_2, 1, x2_frame);

e_x1_est = x1_frame - x1_est;
e_x2_est = x2_frame - x2_est;

% pre-cursor coefficient estimation
L          = length(s_frame);
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


%% Generate/Plot IR for total DAP system (pass impulse through)

s_impulse = zeros(100,1);
s_impulse(1) = 1;

h_channel_1 = filter(b_air_1, a_air_1, s_impulse);
h_channel_2 = filter(b_air_2, a_air_2, s_impulse);

x1_frame = filter(b_air_1, a_air_1, s_impulse);
x2_frame = filter(b_air_2, a_air_2, s_impulse);
h_1_dap = x1_frame - ...
          ((filter(h_pred_1_from_1, 1, x1_frame) + ...
          filter(h_pred_1_from_2, 1, x2_frame)));

h_2_dap = x2_frame - ...
          ((filter(h_pred_2_from_1, 1, x1_frame) + ...
          filter(h_pred_2_from_2, 1, x2_frame)));

h_dap = h_1_dap .* h0(1) + h_2_dap .* h0(2);


Nfft    = 4096;

H_1_dap = fft(h_1_dap, Nfft);
Hm_1_dap  = abs(H_1_dap(1:(Nfft/2)));
Hp_1_dap  = angle(H_1_dap(1:(Nfft/2)));

H_2_dap = fft(h_2_dap, Nfft);
Hm_2_dap  = abs(H_2_dap(1:(Nfft/2)));
Hp_2_dap  = angle(H_2_dap(1:(Nfft/2)));

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

figure()
subplot(5,1,1)
stem((0:(length(h_channel_1)-1)) .* (1/fs), h_channel_1);
xlabel('Time [sec]')
title("Channel 1 Impulse Response")
subplot(5,1,1)
stem((0:(length(h_channel_2)-1)) .* (1/fs), h_channel_2);
xlabel('Time [sec]')
title("Channel 2 Impulse Response")
subplot(5,1,3)
stem((0:(length(h_1_dap)-1)) .* (1/fs), h_1_dap);
xlabel('Time [sec]')
title("h 1 dap Equalized Impulse Response")
subplot(5,1,4)
stem((0:(length(h_2_dap)-1)) .* (1/fs), h_2_dap);
xlabel('Time [sec]')
title("h 2 dap Equalized Impulse Response")
subplot(5,1,5)
stem((0:(length(h_dap)-1)) .* (1/fs), h_dap);
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
title("Equalized Impulse Response")


figure()
subplot(2,1,1)
plot(f_H_air ./ 1000, 20*log10(Hm_channel_1));
title('Channel 1 magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(f_H_air ./ 1000, Hp_channel_1);
title('Channel 1 phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

figure()
subplot(2,1,1)
plot(f_H_air ./ 1000, 20*log10(Hm_channel_2));
title('Channel 2 magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(f_H_air ./ 1000, Hp_channel_1);
title('Channel 2 phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

figure()
subplot(1,2,1)
zplane(b_air_1', a_air_1')
title('Channel 1')
subplot(1,2,2)
zplane(b_air_2', a_air_2')
title('Channel 2')

figure()
subplot(2,1,1)
plot(f_H_air ./ 1000, 20*log10(Hm_1_dap));
title('Hm 1 dap Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(f_H_air ./ 1000, Hp_1_dap);
title('Hp 1 dap Equalized phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

figure()
subplot(3,2,1)
plot(f_H_air ./ 1000, 20*log10(Hm_1_dap));
title('Hm 1 dap Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(3,2,2)
plot(f_H_air ./ 1000, unwrap(Hp_1_dap));
title('Hp 1 dap Equalized phase response')
xlabel('Frequency [kHz]')
ylabel('rad')
subplot(3,2,3)
plot(f_H_air ./ 1000, 20*log10(Hm_2_dap));
title('Hm 2 dap Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(3,2,4)
plot(f_H_air ./ 1000, unwrap(Hp_2_dap));
title('Hp 2 dap Equalized phase response')
xlabel('Frequency [kHz]')
ylabel('rad')
subplot(3,2,5)
plot(f_H_air ./ 1000, 20*log10(Hm_dap));
title('Hm dap Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(3,2,6)
plot(f_H_air ./ 1000, unwrap(Hp_dap));
title('Hp dap Equalized phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

% %% Energy Decay Curve Comparison
% 
% % energy_air_1 = h_channel_1' * h_channel_1;
% % energy_dap   = h_dap' * h_dap;
% % edc_air_1    = zeros(length(h_channel_1), 1);
% % edc_dap      = zeros(length(h_channel_1), 1);
% % for n = 1:length(h_channel_1)
% % 
% %     edc_air_1(n) = energy_air_1;
% %     edc_dap(n)   = energy_dap;
% % 
% %     energy_air_1 = energy_air_1 - h_channel_1(n)^2;
% %     energy_dap   = energy_dap - h_dap(n)^2;
% % end
% 
% figure()
% plot((0:(length(h_channel_1)-1)) .* (1/fs), 20*log10(edc_air_1));
% hold on;
% plot((0:(length(h_channel_1)-1)) .* (1/fs), 20*log10(edc_dap .* (max(edc_air_1) / max(edc_dap))));
% xlabel('Time [sec]')
% ylabel('Instantaneous power [dB]')
% legend('Reverb Energy', 'Equalized Reverb Energy')
% title('Energy Decay Curve')
% % figure()
% % plot((0:(length(h_air_1)-1)) .* (1/fs), 20*log10(energy_air_1));
% % hold on;
% % plot((0:(length(h_dap)-1)) .* (1/fs), 20*log10(energy_dap .* (max(energy_air_1) / max(energy_dap))));
% % xlabel('Time [sec]')
% % ylabel('Energy [dB]')
% % legend('Reverb Energy', 'Equalized Reverb Energy')
% % title('Energy Decay Curve')

