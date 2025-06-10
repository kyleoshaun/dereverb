function [Y, refstim] = compute_mic_signals(source_data, noise_data, interf_data, fs)
% Inputs
% - source_data.signal:   Clean Speech Signal
% - source_data.rir_data: M-channel RIR data for source (cols = RIRs, M cols) 
% - source_data.memo:     Descriptor string
% - source_data.rir_memo: Short Descriptor for room
% - source_data.rir_desc: Long Descriptor for room
% - source_data.stimdb:   dB SPL at ear drum
% 
% - noise_data.enable:    Enable Noise
% - noise_data.signals:   M-channel noise signal data (cols = signals, M cols) -- cols length must match clean speech
% - noise_data.SNR_dB:   Signal-to-noise ratio (dB)
% - noise_data.memo:     Descriptor string
%
% - interf_data.enable:   Enable interference signal (secondary talker)
% - interf_data.signal:   Secondary clean speech signal (length must match clean speech)
% - interf_data.rir_data: M-channel RIR data for interfering talker 
% - interf_data.SIR_dB:      Signal-to-interference ratio (dB)
% - interf_data.memo:     Descriptor string
%
% - fs: Sample rate [Hz]
%
% Outputs
% - Y: Microphone Signals (S*G_s + I*G_i + N) in Pa
% - refstim: Clean source signal calibrated to stimdb dB SPL
%
% Note:
% - Channel Indexing convention: ch1 = Left1, ch2 = right1, ch3=left2, ch4=right2, ...

%% Parameters

% Extract Source Data
s_source       = source_data.signal;
G_source       = source_data.rir_data;
[L_channel, M] = size(G_source);

%% Calibrate speech level (dbstim dB SPL)

% Calibrate clean reference signal level
refstim = source_data.signal / rms(source_data.signal) * 20e-6*10^(source_data.stimdb/20);

% Compute calibration gain for reverberant signal using source * early reflections
ER_length   = round((50 / 1000) * fs); % first 50 msec are selected as "early reflections"
ER_length   = min(ER_length, length(G_source(:,1)));
G_source_ER = G_source(1:ER_length, :);
Y_ER = zeros(length(s_source), M);
for ch = 1:M
    Y_ER(:,ch) = filter(G_source_ER(:,ch), 1, s_source);
end
max_rms_y_ER   = max(rms(Y_ER));
g_calib_source = 20e-6*10^(source_data.stimdb/20) / max_rms_y_ER;
Y_ER  = Y_ER  .* g_calib_source;
max_rms_y_ER = max(rms(Y_ER));

% Apply calibration gain to clean signal before using it to compute microphone signals
s_source = s_source .* g_calib_source;

%% Compute reverberant signals
Y = zeros(length(s_source), M);
for ch = 1:M
    Y(:,ch) = filter(G_source(:,ch), 1, s_source);
end

%% Compute noise/interference levels

if noise_data.enable
    % Extract and validate Noise data
    noise_signals = noise_data.signals;
    [~, M_n]      = size(noise_signals);
    SNR_dB        = noise_data.SNR_dB;

    if M_n ~= M
        error('# of channels for source RIRs and noise data must match')
    end

    if length(noise_signals(:,1)) < length(s_source)
        error('Not enough noise data')
    end

    % Match noise length to source length
    noise_signals = noise_signals(1:length(s_source), :);

    % Set up SNR based on source * early reflections
    max_noise_rms = max(rms(noise_signals));
    snr_i          = max_rms_y_ER / max_noise_rms;
    snr_f          = 10 ^ (SNR_dB / 20);
    noise_gain     = 1 / (snr_f/snr_i);
    noise_signals  = noise_signals .* noise_gain;
end

if interf_data.enable
    % Extract and validate interference data
    interf_signal = interf_data.signal;
    G_interf      = interf_data.rir_data;
    [~, M_i]      = size(G_interf);
    sir_db        = interf_data.SIR_dB;

    if M_i ~= M
        error('# of channels for source and interference RIR data must match')
    end

    if length(interf_signal) ~= length(s_source)
        error('Not enough data for interference signal')
    end

    % Match interference length to source length
    interf_signal = interf_signal(1:length(s_source), :);

    % Compute Reverberant interference
    interf_signals = zeros(length(s_source), M);
    for ch = 1:M
        interf_signals(:,ch) = filter(G_interf(:,ch), 1, interf_signal);
    end

    % Set up SIR based on source * early reflections
    rms_interf     = max(rms(interf_signals));
    rms_y          = max(rms(Y));
    sir_i          = rms_y / rms_interf;
    sir_f          = 10 ^ (sir_db / 20);
    interf_gain    = 1 / (sir_f/sir_i);
    interf_signals  = interf_signals .* interf_gain;
end

%% Compute Microphone Signals

if noise_data.enable
    % Add noise
    Y = Y + noise_signals;
end

if interf_data.enable
    % Add Interference
    Y = Y + interf_signals;
end


end