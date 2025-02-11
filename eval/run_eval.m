function run_eval(M, enable_source_whitening)
% Inputs
% - M: Number of Microphones
% - enable_source_whitening: Set to true to enable stage 1 of delay-and-predict (source whitening)

%% Notes
%
% - All Plots/Metrics describing the "equalized response" or "energy decay"
%   should be taken with a grain of salt since overwhitening of the speech signal
%   due to inaccuracies in the source whitening stage will appear as un-equalized
%   reverb when really it is not
%
% TODO: 
% - Figure out why the time alignment isnt working (Add sample delay to 
%   2-channel setup with L_channel = 10 to see).
%   --> Update: Definitely seems to be broken, also in old scripts --
%               review delay filter updates in original paper
%   --> Update: Even with manual time alignment of the RIRs, 
%               it doesnt seem to work, is this really to do with time alignment at all?
% - Save log to file
% - Save figures to file
%
% Questions
% - Why does MC-LP break down over a certain T60? Does this maybe correspond
%   to channels becoming non-minimimum phase (Note however that just because
%   the individual channels are non-minimum phase doenst mean the MINT filter would
%   include non-minimum phase EQs -- should check this). If MINT includes
%   non-minimum phase EQs, would benefit from covariance method and/or delay.
% - Why does reverb get worse in some cases when running on real speech
%   (not just overwhiten), it seems weird because when S1 is disabled 
%   it works well so whats the diff? Is this just the conditioning of R_mc 
%   due to how white the source is? Does this get better by regularizing R_mc?
%   What about S1 on clean speech?
%
%

addpath('../../samples/')
addpath('../../RIR_Databases/RIR_1_4_BinauralDatabase/')
addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/')
addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/')
addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/courtyard/')
addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/cafeteria/')
addpath('../../RIR_Databases/MYRiAD_V2/tools/MATLAB')
addpath('../../dereverb/LPC/')
addpath('../utilities/matlab')

% Shuffle and save seed for reproducibility
rng('shuffle')
seed = abs(round(randn(1,1) * 10000));
rng(seed);
rng(724); % Reproduce specific seed

%% Parameters

num_cpus = 10;

adjust_gain_ambiguity   = true;
SNR_dB                  = 300;

% Choose Source signal
source = "SA1.WAV";
%[s,fs] = audioread("assgn1.wav");
%[s,fs] = audioread("SA1_phoneme.wav");

[s,fs] = audioread(source);

% Loop to synthetically increase length
% NOTE: MC-LPC stage performance is highly influenced by length of signal
%       Even if stage 1 is bypassed (For very long signals stage 1 bypassed
%       Results in nearly zero forcing)
num_loops_s = round(8.8*fs / length(s));
s = repmat(s, num_loops_s, 1);

if enable_source_whitening == false
    % Set s to white noise (rather than speech) to skip source whitening and
    % debug MC-LPC alone
    s = randn(8.8*fs,1); 
    % Note: The longer the channel IR, the more data is needed,
    % Length 5000000 will give near zero forcing for L_channel = 1000 (This is 312 sec at fs=16kHz though)
    source = "white noise";
end

% Generate randomly shaped noise for speech
% s = randn(140800,1); 
% h_speech = randn(4,1);
% s = filter(h_speech, 1, s);


% FFT/PSD params for plots
Nfft = 4096;
Nfft_welch  = 2^(nextpow2(length(s)));
Nfft_welch  = min(Nfft_welch, Nfft);
Nwin_welch  = Nfft_welch / 2;
Nover_welch = Nfft_welch / 4;

[Sm, w_welch] = pwelch(s, Nwin_welch, Nover_welch, Nfft_welch);
freqs_welch         = w_welch .* (fs / (2*pi));

%% Channel

channel_memo = "synthetic (exponentially decaying white noise)";

% Exponential decay curve
T60 = 250 / 1000;
N60 = T60 * fs;
L_channel = round(N60*1.5);
%tau = L_channel/4;
tau = N60 / log(10^1.5); % -ln(10^(-30/20)) = ln(10^1.5)
exp_decay = exp(-1 .* (1:L_channel)' ./ tau);

% Option 1: Generate random RIR
b_rir_1 = randn(L_channel,1);
b_rir_2 = randn(L_channel,1);
b_rir_3 = randn(L_channel,1);
b_rir_4 = randn(L_channel,1);
b_rir_5 = randn(L_channel,1);
b_rir_6 = randn(L_channel,1);
b_rir_7 = randn(L_channel,1);
b_rir_8 = randn(L_channel,1);

a_rir_1 = 1;
a_rir_2 = 1;
a_rir_3 = 1;
a_rir_4 = 1;
a_rir_5 = 1;
a_rir_6 = 1;
a_rir_7 = 1;
a_rir_8 = 1;

% Apply synthetic exponential decay to RIRs
b_rir_1 = b_rir_1(1:L_channel) .* exp_decay;
b_rir_2 = b_rir_2(1:L_channel) .* exp_decay;
b_rir_3 = b_rir_3(1:L_channel) .* exp_decay;
b_rir_4 = b_rir_4(1:L_channel) .* exp_decay;
b_rir_5 = b_rir_5(1:L_channel) .* exp_decay;
b_rir_6 = b_rir_6(1:L_channel) .* exp_decay;
b_rir_7 = b_rir_7(1:L_channel) .* exp_decay;
b_rir_8 = b_rir_8(1:L_channel) .* exp_decay;

% Synthetic 2-sided RIR
%b_rir_1 = [flip(b_rir_1) ; b_rir_1];
%b_rir_2 = [flip(b_rir_2) ; b_rir_2];

% % Add Synthetic time delay
delay_1 = 0;
delay_2 = 0;
delay_3 = 0;
delay_4 = 0;
delay_5 = 0;
delay_6 = 0;
delay_7 = 0;
delay_8 = 0;
max_delay = max([delay_1 delay_2 delay_3 delay_4 delay_5 delay_6 delay_7 delay_8]);
b_rir_1 = [zeros(delay_1,1) ; b_rir_1 ; zeros((max_delay - delay_1),1)];
b_rir_2 = [zeros(delay_2,1) ; b_rir_2 ; zeros((max_delay - delay_2),1)];
b_rir_3 = [zeros(delay_3,1) ; b_rir_3 ; zeros((max_delay - delay_3),1)];
b_rir_4 = [zeros(delay_4,1) ; b_rir_4 ; zeros((max_delay - delay_4),1)];
b_rir_5 = [zeros(delay_5,1) ; b_rir_5 ; zeros((max_delay - delay_5),1)];
b_rir_6 = [zeros(delay_6,1) ; b_rir_6 ; zeros((max_delay - delay_6),1)];
b_rir_7 = [zeros(delay_7,1) ; b_rir_7 ; zeros((max_delay - delay_7),1)];
b_rir_8 = [zeros(delay_8,1) ; b_rir_8 ; zeros((max_delay - delay_8),1)];

rir_memo = "synthRIR";

real_rir = 0;
RIR_database = 1; % 0 = HRIR, 1 = MYRiAD
if real_rir
rir_memo = "realRIR";
channel_memo = "Measured RIR";
% Option 2: Real RIR

if RIR_database == 0 % HRIR Database

    if M > 6
        error("HRIR RIR Database only supports up to 6 microphones")
    end

    % Channels:
    % - 1: Left Front
    % - 2: Right Front
    % - 3: Left Middle
    % - 4: Right Middle
    % - 5: Left Rear
    % - 6: Front Rear
    head_orientation = 1;
    speaker_loc      = 'A';
    data_set         = 'bte';
    HRIR_data = loadHRIR('Courtyard', head_orientation, speaker_loc, data_set);
    if HRIR_data.fs ~= fs
        [resample_p, resample_q] = rat(fs / HRIR_data.fs);
        HRIR_data.data = resample(HRIR_data.data, resample_p, resample_q);
    end
    
    b_rir_1 = HRIR_data.data(:,1);
    b_rir_2 = HRIR_data.data(:,2);
    b_rir_3 = HRIR_data.data(:,3);
    b_rir_4 = HRIR_data.data(:,4);
    b_rir_5 = HRIR_data.data(:,5);
    b_rir_6 = HRIR_data.data(:,6);
    b_rir_7 = b_rir_5;
    b_rir_8 = b_rir_6;

    edc_1 = EDC(b_rir_1);
    edc_2 = EDC(b_rir_2);
    edc_3 = EDC(b_rir_3);
    edc_4 = EDC(b_rir_4);
    edc_5 = EDC(b_rir_5);
    edc_6 = EDC(b_rir_6);
    
    figure()
    subplot(1,2,1)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_1)
    hold on;
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_2)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_3)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_4)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_5)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_6)
    xlabel('Time [sec]')
    legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5', 'Channel 6')
    title(sprintf('Real Measured RIRs (HRIR %s)', room))
    
    subplot(1,2,2)
    plot((0:(length(edc_1)-1)) .* (1/fs), 10*log10(edc_1));
    hold on;
    plot((0:(length(edc_2)-1)) .* (1/fs), 10*log10(edc_2));
    plot((0:(length(edc_3)-1)) .* (1/fs), 10*log10(edc_3));
    plot((0:(length(edc_4)-1)) .* (1/fs), 10*log10(edc_4));
    plot((0:(length(edc_5)-1)) .* (1/fs), 10*log10(edc_5));
    plot((0:(length(edc_6)-1)) .* (1/fs), 10*log10(edc_6));
    grid on;
    ylim([-65 6])
    xlabel('Time [sec]')
    ylabel('dB')
    legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5', 'Channel 6')
    title(sprintf('Energy Decay Curve (HRIR %s)', room))

else

    if M > 4
        error("MYRiAD RIR Database only supports up to 4 microphones")
    end


    MYRiAD = my_load_audio_data;
    MYRiAD_audio = MYRiAD.data{1};
    fs_MYRiAD = 44100;

    if fs_MYRiAD ~= fs
        [resample_p, resample_q] = rat(fs / fs_MYRiAD);
        MYRiAD_audio = resample(MYRiAD_audio, resample_p, resample_q);
    end
    
    b_rir_1 = MYRiAD_audio(:,1);
    b_rir_2 = MYRiAD_audio(:,2);
    b_rir_3 = MYRiAD_audio(:,3);
    b_rir_4 = MYRiAD_audio(:,4);
    b_rir_5 = b_rir_1;
    b_rir_6 = b_rir_2;
    b_rir_7 = b_rir_3;
    b_rir_8 = b_rir_4;
    
    edc_1 = EDC(b_rir_1);
    edc_2 = EDC(b_rir_2);
    edc_3 = EDC(b_rir_3);
    edc_4 = EDC(b_rir_4);
    
    figure()
    subplot(1,2,1)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_1)
    hold on;
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_2)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_3)
    plot((0:(length(b_rir_1)-1)) .* (1/fs), b_rir_4)
    xlabel('Time [sec]')
    legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4')
    title(sprintf('Real Measured RIRs (MYRiAD SAL)'))
    
    subplot(1,2,2)
    plot((0:(length(edc_1)-1)) .* (1/fs), 10*log10(edc_1));
    hold on;
    plot((0:(length(edc_2)-1)) .* (1/fs), 10*log10(edc_2));
    plot((0:(length(edc_3)-1)) .* (1/fs), 10*log10(edc_3));
    plot((0:(length(edc_4)-1)) .* (1/fs), 10*log10(edc_4));
    grid on;
    ylim([-65 6])
    xlabel('Time [sec]')
    ylabel('dB')
    legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4')
    title(sprintf('Energy Decay Curve (MYRiAD SAL)'))

end

% Remove measurement noise manually (dont want to convolve this, it isnt really part of the RIR)
meas_noise_db_thresh = 15; % -15 dB considered measurement noise
[b_rir_1, tof_1_est] = remove_air_meas_noise(b_rir_1, meas_noise_db_thresh);
[b_rir_2, tof_2_est] = remove_air_meas_noise(b_rir_2, meas_noise_db_thresh);
[b_rir_3, tof_3_est] = remove_air_meas_noise(b_rir_3, meas_noise_db_thresh);
[b_rir_4, tof_4_est] = remove_air_meas_noise(b_rir_4, meas_noise_db_thresh);
[b_rir_5, tof_5_est] = remove_air_meas_noise(b_rir_5, meas_noise_db_thresh);
[b_rir_6, tof_6_est] = remove_air_meas_noise(b_rir_6, meas_noise_db_thresh);
[b_rir_7, tof_7_est] = remove_air_meas_noise(b_rir_7, meas_noise_db_thresh);
[b_rir_8, tof_8_est] = remove_air_meas_noise(b_rir_8, meas_noise_db_thresh);

% Compensate Delays
tofs_est   = [tof_1_est tof_2_est tof_3_est tof_4_est tof_5_est tof_6_est tof_7_est tof_8_est];
max_tof    = max(tofs_est);
delays_tof = max_tof - tofs_est;
max_delay  = max_tof - min(tofs_est);
b_rir_1    = [zeros(delays_tof(1),1) ; b_rir_1 ; zeros((max_delay - delays_tof(1)),1)];
b_rir_2    = [zeros(delays_tof(2),1) ; b_rir_2 ; zeros((max_delay - delays_tof(2)),1)];
b_rir_3    = [zeros(delays_tof(3),1) ; b_rir_3 ; zeros((max_delay - delays_tof(3)),1)];
b_rir_4    = [zeros(delays_tof(4),1) ; b_rir_4 ; zeros((max_delay - delays_tof(4)),1)];
b_rir_5    = [zeros(delays_tof(5),1) ; b_rir_5 ; zeros((max_delay - delays_tof(5)),1)];
b_rir_6    = [zeros(delays_tof(6),1) ; b_rir_6 ; zeros((max_delay - delays_tof(6)),1)];
b_rir_7    = [zeros(delays_tof(7),1) ; b_rir_7 ; zeros((max_delay - delays_tof(7)),1)];
b_rir_8    = [zeros(delays_tof(8),1) ; b_rir_8 ; zeros((max_delay - delays_tof(8)),1)];


% Truncate to shorten
% trunc_length = length(b_rir_1);
% tof_bulk = 1;
% b_rir_1 = b_rir_1(tof_bulk:trunc_length);
% b_rir_2 = b_rir_2(tof_bulk:trunc_length);
% b_rir_3 = b_rir_3(tof_bulk:trunc_length);
% b_rir_4 = b_rir_4(tof_bulk:trunc_length);
% b_rir_5 = b_rir_5(tof_bulk:trunc_length);
% b_rir_6 = b_rir_6(tof_bulk:trunc_length);

% figure()
% subplot(6,1,1)
% plot(b_rir_1)
% subplot(6,1,2)
% plot(b_rir_2)
% subplot(6,1,3)
% plot(b_rir_3)
% subplot(6,1,4)
% plot(b_rir_4)
% subplot(6,1,5)
% plot(b_rir_5)
% subplot(6,1,6)
% plot(b_rir_6)
% sgtitle('Real Measured RIRs after conditioning')

end

auto_time_align = 0;
if auto_time_align

fprintf("Manually Time-Aligning RIRs...\n")

[r_21, lags_21] = xcorr(b_rir_2, b_rir_1);
[r_31, lags_31] = xcorr(b_rir_3, b_rir_1);
[r_41, lags_41] = xcorr(b_rir_4, b_rir_1);
[r_51, lags_51] = xcorr(b_rir_5, b_rir_1);
[r_61, lags_61] = xcorr(b_rir_6, b_rir_1);
[r_71, lags_71] = xcorr(b_rir_7, b_rir_1);
[r_81, lags_81] = xcorr(b_rir_8, b_rir_1);

tofs_est = [ ...
0 ...
lags_21(find(abs(r_21) == max(abs(r_21)))) ...
lags_31(find(abs(r_31) == max(abs(r_31)))) ...
lags_41(find(abs(r_41) == max(abs(r_41)))) ...
lags_51(find(abs(r_51) == max(abs(r_51)))) ...
lags_61(find(abs(r_61) == max(abs(r_61)))) ...
lags_71(find(abs(r_71) == max(abs(r_71)))) ...
lags_81(find(abs(r_81) == max(abs(r_81)))) ...
];

fprintf("  - Estimated Time of Flights Before Alignment = [")
for ch = 1:M
    fprintf("%.0f ", tofs_est(ch));
end
fprintf("]\n");

total_avg_time_removed = 0;
total_avg_time_added   = 0;
%while (max(abs(tofs)) ~= 0) % Perfectly aligned
while (sum(abs(tofs_est) ~= 0) > 2) % Most channels are perfectly aligned
    
    % Manual time alignment via addition of delay
    max_tof = max(tofs_est);
    delays_tof  = max_tof - tofs_est;
    b_rir_1 = [zeros(delays_tof(1),1) ; b_rir_1(1:(end-delays_tof(1)))];
    b_rir_2 = [zeros(delays_tof(2),1) ; b_rir_2(1:(end-delays_tof(2)))];
    b_rir_3 = [zeros(delays_tof(3),1) ; b_rir_3(1:(end-delays_tof(3)))];
    b_rir_4 = [zeros(delays_tof(4),1) ; b_rir_4(1:(end-delays_tof(4)))];
    b_rir_5 = [zeros(delays_tof(5),1) ; b_rir_5(1:(end-delays_tof(5)))];
    b_rir_6 = [zeros(delays_tof(6),1) ; b_rir_6(1:(end-delays_tof(6)))];
    b_rir_7 = [zeros(delays_tof(7),1) ; b_rir_7(1:(end-delays_tof(7)))];
    b_rir_8 = [zeros(delays_tof(8),1) ; b_rir_8(1:(end-delays_tof(8)))];
    total_avg_time_added = total_avg_time_added + mean(delays_tof);
    
    % Manual time alignment via removal of delay
    % min_tof_est  = min(tofs_est);
    % rel_tofs_est = tofs_est - min_tof_est; 
    % b_rir_1 = [b_rir_1((rel_tofs_est(1)+1):end) ; zeros(rel_tofs_est(1),1)];
    % b_rir_2 = [b_rir_2((rel_tofs_est(2)+1):end) ; zeros(rel_tofs_est(2),1)];
    % b_rir_3 = [b_rir_3((rel_tofs_est(3)+1):end) ; zeros(rel_tofs_est(3),1)];
    % b_rir_4 = [b_rir_4((rel_tofs_est(4)+1):end) ; zeros(rel_tofs_est(4),1)];
    % b_rir_5 = [b_rir_5((rel_tofs_est(5)+1):end) ; zeros(rel_tofs_est(5),1)];
    % b_rir_6 = [b_rir_6((rel_tofs_est(6)+1):end) ; zeros(rel_tofs_est(6),1)];
    % b_rir_7 = [b_rir_7((rel_tofs_est(7)+1):end) ; zeros(rel_tofs_est(7),1)];
    % b_rir_8 = [b_rir_8((rel_tofs_est(8)+1):end) ; zeros(rel_tofs_est(8),1)];
    % total_avg_time_removed = total_avg_time_removed + mean(rel_tofs_est);
    
    
    % Re-estimate TOFs to show results
    [r_21, lags_21] = xcorr(b_rir_2, b_rir_1);
    [r_31, lags_31] = xcorr(b_rir_3, b_rir_1);
    [r_41, lags_41] = xcorr(b_rir_4, b_rir_1);
    [r_51, lags_51] = xcorr(b_rir_5, b_rir_1);
    [r_61, lags_61] = xcorr(b_rir_6, b_rir_1);
    [r_71, lags_71] = xcorr(b_rir_7, b_rir_1);
    [r_81, lags_81] = xcorr(b_rir_8, b_rir_1);
    
    tofs_est = [ ...
    0 ...
    lags_21(find(abs(r_21) == max(abs(r_21)))) ...
    lags_31(find(abs(r_31) == max(abs(r_31)))) ...
    lags_41(find(abs(r_41) == max(abs(r_41)))) ...
    lags_51(find(abs(r_51) == max(abs(r_51)))) ...
    lags_61(find(abs(r_61) == max(abs(r_61)))) ...
    lags_71(find(abs(r_71) == max(abs(r_71)))) ...
    lags_81(find(abs(r_81) == max(abs(r_81)))) ...
    ];
    
    fprintf("  - Estimated Time of Flights After Alignment = [")
    for ch = 1:M
        fprintf("%.0f ", tofs_est(ch));
    end
    fprintf("]\n");

    %break;
end

fprintf("Done (Total avg time removed = %.2f Total avg delay added = %.2f\n", total_avg_time_removed, total_avg_time_added);



end


%% Compute Microphone Signals

% Compute reverberant signals
y1 = filter(b_rir_1, a_rir_1, s);
y2 = filter(b_rir_2, a_rir_2, s);
y3 = filter(b_rir_3, a_rir_3, s);
y4 = filter(b_rir_4, a_rir_4, s);
y5 = filter(b_rir_5, a_rir_5, s);
y6 = filter(b_rir_6, a_rir_6, s);
y7 = filter(b_rir_7, a_rir_7, s);
y8 = filter(b_rir_8, a_rir_8, s);

% Add measurement noise
% NOTE: This is essentially regularizing (biasing) R_mc so dont make it unrealistically high
noise_mag = rms(s) .* 10 ^ (-1*SNR_dB / 20);
y1 = y1 + noise_mag .* randn(length(y1), 1);
y2 = y2 + noise_mag .* randn(length(y2), 1);
y3 = y3 + noise_mag .* randn(length(y3), 1);
y4 = y4 + noise_mag .* randn(length(y4), 1);
y5 = y5 + noise_mag .* randn(length(y5), 1);
y6 = y6 + noise_mag .* randn(length(y6), 1);
y7 = y7 + noise_mag .* randn(length(y7), 1);
y8 = y8 + noise_mag .* randn(length(y8), 1);

%% Interfering talker

% [sv,fs] = audioread("SA2.wav");
% 
% num_loops_s = round(10*fs / length(sv));
% sv = repmat(sv, num_loops_s, 1);
% 
% % Option 2: Real RIR
% % Channels:
% % - 1: Left Front
% % - 2: Right Front
% % - 3: Left Middle
% % - 4: Right Middle
% % - 5: Left Rear
% % - 6: Front Rear
% head_orientation = 2;
% speaker_loc      = 'C';
% data_set         = 'bte';
% HRIR_data = loadHRIR('office_II', head_orientation, speaker_loc, data_set);
% if HRIR_data.fs ~= fs
%     [resample_p, resample_q] = rat(fs / HRIR_data.fs);
%     HRIR_data.data = resample(HRIR_data.data, resample_p, resample_q);
% end
% 
% b_rir_v_1 = HRIR_data.data(:,1);
% b_rir_v_2 = HRIR_data.data(:,2);
% b_rir_v_3 = HRIR_data.data(:,3);
% b_rir_v_4 = HRIR_data.data(:,4);
% b_rir_v_5 = HRIR_data.data(:,5);
% b_rir_v_6 = HRIR_data.data(:,6);
% b_rir_v_7 = b_rir_v_5;
% b_rir_v_8 = b_rir_v_6;
% 
% b_rir_v_1 = b_rir_v_1(tof_bulk:trunc_length);
% b_rir_v_2 = b_rir_v_2(tof_bulk:trunc_length);
% b_rir_v_3 = b_rir_v_3(tof_bulk:trunc_length);
% b_rir_v_4 = b_rir_v_4(tof_bulk:trunc_length);
% b_rir_v_5 = b_rir_v_5(tof_bulk:trunc_length);
% b_rir_v_6 = b_rir_v_6(tof_bulk:trunc_length);
% b_rir_v_7 = b_rir_v_7(tof_bulk:trunc_length);
% b_rir_v_8 = b_rir_v_8(tof_bulk:trunc_length);
% 
% % Compute reverberant signals
% v1 = filter(b_rir_v_1, 1, sv);
% v2 = filter(b_rir_v_2, 1, sv);
% v3 = filter(b_rir_v_3, 1, sv);
% v4 = filter(b_rir_v_4, 1, sv);
% v5 = filter(b_rir_v_5, 1, sv);
% v6 = filter(b_rir_v_6, 1, sv);
% v7 = filter(b_rir_v_7, 1, sv);
% v8 = filter(b_rir_v_8, 1, sv);
% 
% % Add measurement noise
% % NOTE: This is essentially regularizing (biasing) R_mc so dont make it unrealistically high
% noise_mag = rms(sv) .* 10 ^ (-1*SNR_dB / 20);
% v1 = v1 + noise_mag .* randn(length(v1), 1);
% v2 = v2 + noise_mag .* randn(length(v2), 1);
% v3 = v3 + noise_mag .* randn(length(v3), 1);
% v4 = v4 + noise_mag .* randn(length(v4), 1);
% v5 = v5 + noise_mag .* randn(length(v5), 1);
% v6 = v6 + noise_mag .* randn(length(v6), 1);
% v7 = v7 + noise_mag .* randn(length(v7), 1);
% v8 = v8 + noise_mag .* randn(length(v8), 1);
% 
% % y1 = [y1 ; v1];
% % y2 = [y2 ; v2];
% % y3 = [y3 ; v3];
% % y4 = [y4 ; v4];
% % y5 = [y5 ; v5];
% % y6 = [y6 ; v6];
% % y7 = [y7 ; v7];
% % y8 = [y8 ; v8];
% % 
% % s = [s ; sv];

%% Run Delay-and-Predict

% % TODO: Test why this isnt working
% b_rir_1 = [0 ; b_rir_1];
% b_rir_2 = [b_rir_2 ; 0];

%b_rir_1 = [b_rir_1(2:end) ; 0];

switch M
    case 2
        Y = [y1 y2];
    case 3
        Y = [y1 y2 y3];
    case 4
        Y = [y1 y2 y3 y4];
    case 5
        Y = [y1 y2 y3 y4 y5];
    case 6
        Y = [y1 y2 y3 y4 y5 y6];
    case 7
        Y = [y1 y2 y3 y4 y5 y6 y7];
    case 8
        Y = [y1 y2 y3 y4 y5 y6 y7 y8];
end

% Prediction Orders
L_channel = length(b_rir_1);
p2 = round(L_channel * 1.25 / (M-1)); % Stage 2 MC-LPC order
p1 = 2*(L_channel+p2);%round(p2 + L_channel * 1.25); % Stage 1 Source Whitening order

% Create Results Folder
s1_memo = "1s1";
if enable_source_whitening == false
    s1_memo = "0s1";
end
results_dir = sprintf('Results/dap_%s_%.0fM_%s_%.0fL_%.0fp1_%.0fp2_%.0fkFs_%s', rir_memo, M, s1_memo, L_channel, p1, p2, fs / 1000, datetime('today'));
if ~exist(results_dir,'Dir')
    mkdir(results_dir)
end

% Start Logging Console if exist(
diary_fullpath = sprintf('%s/console.log', results_dir);
if exist(diary_fullpath, 'file')
    delete(diary_fullpath)
end
diary(diary_fullpath)

fprintf("\nTest Conditions:\n")
fprintf(" - Source Signal = %s\n", source)
fprintf(" - RIR = %s\n", channel_memo)
fprintf(" - L_channel = %.0f\n", L_channel)
fprintf(" - SNR = %.0f\n", SNR_dB)


if num_cpus > 1
    [s_est, dap_equalizer_struct] = delay_and_predict(Y, s, p1, p2, fs, enable_source_whitening, results_dir, num_cpus);
else
    [s_est, dap_equalizer_struct] = delay_and_predict_nopar(Y, s, p1, p2, fs, enable_source_whitening, results_dir);
end

% Adjust for gain/scaling ambiguity
if adjust_gain_ambiguity
    s_est = (rms(s) / rms(s_est)) .* s_est;
end

H      = dap_equalizer_struct.H;
delays = dap_equalizer_struct.delays;
h0     = dap_equalizer_struct.h0;
delay_comp = dap_equalizer_struct.delay_comp;

%% MINT: Ideal solution based on known RIRs

run_MINT = 1;
if run_MINT
    fprintf("Running MINT ... ");
    tau = 0;
    
    switch M
        case 2
            G = [b_rir_1 b_rir_2];
        case 3
            G = [b_rir_1 b_rir_2 b_rir_3];
        case 4
            G = [b_rir_1 b_rir_2 b_rir_3 b_rir_4];
        case 5
            G = [b_rir_1 b_rir_2 b_rir_3 b_rir_4 b_rir_5];
        case 6
            G = [b_rir_1 b_rir_2 b_rir_3 b_rir_4 b_rir_5 b_rir_6];
        case 7
            G = [b_rir_1 b_rir_2 b_rir_3 b_rir_4 b_rir_5 b_rir_6 b_rir_7];
        case 8
            G = [b_rir_1 b_rir_2 b_rir_3 b_rir_4 b_rir_5 b_rir_6 b_rir_7 b_rir_8];
    end
    
    H_mint = MINT(G, tau, results_dir);
    
    save(sprintf('%s/H_mint.mat', results_dir),'H_mint');

    fprintf("Done")

end


%% Plots

%s_delayed  = filter([zeros(delay_comp,1) ; 1], 1, s);

%s_est = apply_dap_equalizer(dap_equalizer_struct, Y);
%s_est = (rms(s) / rms(s_est)) .* s_est;


e_est_dap = s - s_est;

Sm          = pwelch(s,   Nwin_welch, Nover_welch, Nfft_welch);
S_est_dap_m = pwelch(s_est, Nwin_welch, Nover_welch, Nfft_welch);
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

%% Generate/Plot IR for total DAP system (pass impulse through)

s_impulse = zeros(L_channel + p1 + p2,1);
s_impulse(1) = 1;

g_channel_1 = filter(b_rir_1, a_rir_1, s_impulse);
g_channel_2 = filter(b_rir_2, a_rir_2, s_impulse);
g_channel_3 = filter(b_rir_3, a_rir_3, s_impulse);
g_channel_4 = filter(b_rir_4, a_rir_4, s_impulse);
g_channel_5 = filter(b_rir_5, a_rir_5, s_impulse);
g_channel_6 = filter(b_rir_6, a_rir_6, s_impulse);
g_channel_7 = filter(b_rir_7, a_rir_7, s_impulse);
g_channel_8 = filter(b_rir_8, a_rir_8, s_impulse);

switch M
    case 2
        G = [g_channel_1 g_channel_2];
    case 3
        G = [g_channel_1 g_channel_2 g_channel_3];
    case 4
        G = [g_channel_1 g_channel_2 g_channel_3 g_channel_4];
    case 5
        G = [g_channel_1 g_channel_2 g_channel_3 g_channel_4 g_channel_5];
    case 6
        G = [g_channel_1 g_channel_2 g_channel_3 g_channel_4 g_channel_5 g_channel_6];
    case 7
        G = [g_channel_1 g_channel_2 g_channel_3 g_channel_4 g_channel_5 g_channel_6 g_channel_7];
    case 8
        G = [g_channel_1 g_channel_2 g_channel_3 g_channel_4 g_channel_5 g_channel_6 g_channel_7 g_channel_8];
end

eir = apply_dap_equalizer(dap_equalizer_struct, G);

[G_channel_1, freqs_freqz] = freqz(b_rir_1, a_rir_1, Nfft, fs);
Gm_channel_1 = abs(G_channel_1);
Gp_channel_1 = angle(G_channel_1);

G_channel_2  = freqz(b_rir_2, a_rir_2, Nfft);
Gm_channel_2 = abs(G_channel_2);
Gp_channel_2 = angle(G_channel_2);

EIRft = fft(eir, Nfft);
EIRm  = abs(EIRft(1:(Nfft/2)));
EIRp  = angle(EIRft(1:(Nfft/2)));
freqs_fft = (0:(length(EIRm)-1))' .* (fs / Nfft);

figure()
subplot(2,1,1)
stem((0:(length(g_channel_1)-1)) .* (1/fs), g_channel_1);
xlabel('Time [sec]')
title("Channel 1 Impulse Response")
subplot(2,1,2)
stem((0:(length(eir)-1)) .* (1/fs), eir);
xlabel('Time [sec]')
title("Equalized Impulse Response (Stem Plot)")

saveas(gcf, sprintf('%s/EIR_stem.fig', results_dir));

figure()
subplot(2,1,1)
plot((0:(length(g_channel_1)-1)) .* (1/fs), g_channel_1);
xlabel('Time [sec]')
title("Channel 1 Impulse Response")
subplot(2,1,2)
plot((0:(length(eir)-1)) .* (1/fs), eir);
xlabel('Time [sec]')
title("Equalized Impulse Response")

saveas(gcf, sprintf('%s/EIR.fig', results_dir));

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

% figure()
% subplot(1,2,1)
% zplane(b_rir_1', a_rir_1')
% title('Channel 1')
% subplot(1,2,2)
% zplane(b_rir_2', a_rir_2')
% title('Channel 2')

figure()
subplot(2,1,1)
plot(freqs_fft ./ 1000, 20*log10(EIRm));
title('Equalized magnitude Response')
xlabel('Frequency [kHz]')
ylabel('dB')
subplot(2,1,2)
plot(freqs_fft ./ 1000, unwrap(EIRp));
title('Equalized phase response')
xlabel('Frequency [kHz]')
ylabel('rad')

saveas(gcf, sprintf('%s/Equalized_RTF.fig', results_dir));

%% Energy Decay Curve Comparison

edc_rir_1 = EDC(g_channel_1);
edc_dap   = EDC(eir);

figure()
plot((0:(length(g_channel_1)-1)) .* (1/fs), 20*log10(edc_rir_1));
hold on;
plot((0:(length(g_channel_1)-1)) .* (1/fs), 20*log10(edc_dap .* (max(edc_rir_1) / max(edc_dap))));
ylim([-80 60])
%xlim([0 (length(h_channel_1) * (1/fs))])
xlabel('Time [sec]')
ylabel('Energy remaining [dB]')
legend('Reverb Energy', 'Equalized Reverb Energy')
title('Energy Decay Curve')

saveas(gcf, sprintf('%s/EDC.fig', results_dir));

%% Spectrogram Results

% Spectrogram Parameters
N_window  = 256;
N_overlap = N_window / 2;
N_fft     = N_window;

% Compute Spectrograms
window = hamming(N_window);
[S_spec, f_S, t_S] = spectrogram(s, window, N_overlap, N_fft, fs);
[Y_spec, f_Y, t_Y] = spectrogram(y1, window, N_overlap, N_fft, fs);
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
sgtitle('Dereverberation Results')

saveas(gcf, sprintf('%s/Spectrogram.fig', results_dir));

%% Extra plots

enable_extra_plots = 0;
if enable_extra_plots

y1 = filter(b_rir_1, a_rir_1, s_impulse);
y2 = filter(b_rir_2, a_rir_2, s_impulse);

y1 = filter([zeros(delays(1),1) ; 1], 1, y1);
y2 = filter([zeros(delays(2),1) ; 1], 1, y2);

h_pred_1_from_1 = H(:,1);
h_pred_1_from_2 = H(:,2);
h_pred_2_from_1 = H(:,3);
h_pred_2_from_2 = H(:,4);

h_1_dap = y1 - ...
          ((filter(h_pred_1_from_1, 1, y1) + ...
            filter(h_pred_1_from_2, 1, y2)));

h_2_dap = y2 - ...
          ((filter(h_pred_2_from_1, 1, y1) + ...
            filter(h_pred_2_from_2, 1, y2)));
%
% eir = h_1_dap .* h0(1) + h_2_dap .* h0(2);

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

saveas(gcf, sprintf('%s/EIR_2.fig', results_dir));

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

saveas(gcf, sprintf('%s/PredictionFilters.fig', results_dir));

end % Extra plots

%% Test: Apply dereverb to an actual speech signal 

[s,fs] = audioread("SA1.WAV");

% Compute reverberant signals
y1 = filter(b_rir_1, a_rir_1, s);
y2 = filter(b_rir_2, a_rir_2, s);
y3 = filter(b_rir_3, a_rir_3, s);
y4 = filter(b_rir_4, a_rir_4, s);
y5 = filter(b_rir_5, a_rir_5, s);
y6 = filter(b_rir_6, a_rir_6, s);
y7 = filter(b_rir_7, a_rir_7, s);
y8 = filter(b_rir_8, a_rir_8, s);

switch M
    case 2
        Y = [y1 y2];
    case 3
        Y = [y1 y2 y3];
    case 4
        Y = [y1 y2 y3 y4];
    case 5
        Y = [y1 y2 y3 y4 y5];
    case 6
        Y = [y1 y2 y3 y4 y5 y6];
    case 7
        Y = [y1 y2 y3 y4 y5 y6 y7];
    case 8
        Y = [y1 y2 y3 y4 y5 y6 y7 y8];
end

s_est = apply_dap_equalizer(dap_equalizer_struct, Y);

% Spectrogram Parameters
N_window  = 256;
N_overlap = N_window / 2;
N_fft     = N_window;

% Compute Spectrograms
window = hamming(N_window);
[S_spec, f_S, t_S] = spectrogram(s, window, N_overlap, N_fft, fs);
[Y_spec, f_Y, t_Y] = spectrogram(y1, window, N_overlap, N_fft, fs);
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
sgtitle('Dereverberation Results (Reapplied DAP-EQ to SA1.wav)')

saveas(gcf, sprintf('%s/Spectrogram_reappliedToSA1.fig', results_dir));

%% Save results

% Stop logging console
diary off;



end