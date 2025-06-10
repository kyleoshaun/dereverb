restoredefaultpath

close all
clear
clc

%% Add paths

code_folder = "../../";

addpath(sprintf('%s/dereverb/LPC/', code_folder))
addpath(sprintf('%s/dereverb/utilities/matlab', code_folder))
addpath(sprintf('%s/metrics/', code_folder))

addpath(sprintf('%s/samples', code_folder))
addpath(sprintf('%s/RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/', code_folder))
addpath(sprintf('%s/RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/', code_folder))
addpath(sprintf('%s/RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/courtyard/', code_folder))

% Integration of EC with BEZ model
addpath(sprintf("%s/metrics/BSIM2020_EC_BEZ2018a_model/", code_folder))

% HASPI/HASQI
addpath(sprintf("%s/metrics/HASPIv2_HASQIv2_HAAQIv1_Common/", code_folder))

% BSIM
addpath(sprintf('%s/metrics/BSIM_2020/', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/vad/', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/PreProc/', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/gammatonegram/', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/SRMRToolbox-master/libs/auditory/', code_folder))

addpath(genpath(sprintf('%s/metrics/BSIM_2020/Anechoic/', code_folder)))
addpath(sprintf('%s/metrics/BSIM_2020/functions/', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/functions/gammatonefilter/', code_folder))
addpath(sprintf('%s/metrics/BSIM_2020/functions/ltfat/', code_folder))

% Model (BEZ2018a)
addpath(sprintf('%s/metrics/BEZ2018a_model/', code_folder))

% RIR Databases
addpath(sprintf('%s/RIR_Databases', code_folder))

% Other
addpath(sprintf('%s/metrics/NSIM_Papers_ExampleCode 2/NSIMforBEZ2018a_model/', code_folder))
addpath(sprintf('%s/metrics/STMI_Papers_ExampleCode 1/STMIforBEZ2018a_model', code_folder))


%% Shuffle and save seed for reproducibility
rng('shuffle')
seed = abs(round(randn(1,1) * 10000));
rng(seed);
rng(724); % Reproduce specific seed


%% Simulation Parameters

M            = 4; % Number of Microphones
fs           = 16000; % Sample rate (Hz)
stimdb       = 65; % speech level in dB SPL 
enable_noise = false;
SNR_dB       = Inf;
source_length_sec = 10;

% Simulation Flags
enable_debug_plots        = 0;
enable_reduction_of_ERs   = 0;
replace_SAL_with_synth    = 1;
replace_source_with_noise = 0;
s1_enable                 = 1;
s1_on_clean_speech        = 0;

fig_size = [10 10 1300 600];

% Hearing Loss 
HL_freqs = [250 500 1000 2000 4000 6000]; % 6 audiometric frequencies used by HASPI
HL_NH    = [0   0   0    0    0    0]; % Normal Hearing
HL_HI    = [35  35  40   50   60   65]; % High-freq hearing loss (IEC 60118-15 Moderate HL, Moderately Sloping Group)

%% Source Signal

%source_file = 'SA1.WAV';
source_file = 'TMIT_MKLS0.WAV'; % MALE TALKER (all samples concatenated)
% source_file = 'TMIT_MMRP0.WAV'; % MALE TALKER (all samples concatenated)
% source_file = 'TMIT_FAEM0.WAV'; % FEMALE TALKER (TMIT FAEM0, all samples concatenated)
% source_file = 'TMIT_FSAH0.WAV'; % FEMALE TALKER (TMIT FAEM0, all samples concatenated)
% source_file = 'TMIT.WAV'; % Full TMIT Database concatenated

[s_source, Fs_stim] = audioread(source_file);
L_source            = length(s_source);
source_memo         = source_file(1:(end-4)); % remove ".WAV"

if Fs_stim ~= fs
    s_source = resample(s_source, fs, Fs_stim);
end

% Loop to synthetically increase length
num_loops_s = ceil(source_length_sec * fs / length(s_source));
s_source = repmat(s_source, num_loops_s, 1);
s_source = s_source(1:round(source_length_sec * fs));

%% Optional: Substitute white noise for source

if replace_source_with_noise

source_duration_pre_loop = 10;
s_source = randn(source_duration_pre_loop*fs, 1);
num_loops_s = ceil(source_length_sec*fs / length(s_source));
s_source = repmat(s_source, num_loops_s, 1);
source_memo = sprintf("white noise (%d sec looped to %d sec)", source_duration_pre_loop, source_length_sec);

%s_source = randn(length(source_length_sec), 1);
%source_memo = "white noise";

% Apply a custom shape to spectrum
%h_source_shape = 0.65 * conv(conv([1 -0.65 0], [1 0 0.45]), [1 0.5]);
%h_source_shape = 0.65 * conv(conv([1 -0.93 0], [1 0 0.875]), [1 0.5]);
h_source_shape = 0.65 * conv(conv([1 -0.995 0], [1 0 0.99]), [1 0.5]);
s_source = filter(h_source_shape, 1, s_source);
source_memo = sprintf("%s shaped with filter", source_memo);

end

%% Load real RIR data

rir_database = load('rir_databases_4micCombined_1noise.mat');
rir_database = rir_database.rir_databases_4micCombined_1noise;

if rir_database.fs ~= fs
    error("RIR Database sample rate mismatch");
end

room_num       = 4; % SAL (T60 = 2.1 sec)
G_source       = rir_database.rir_data.rir_list_0deg{room_num};
tof_list       = rir_database.rir_data.tof_est_0deg{room_num};
room_memo      = "SALtrunc";
rir_desc       = sprintf("MYRiAD SAL Measured RIR (T60 = 2100 msec, Truncated Exponentially)");

% Manually align RIRs
figure()
for ch=1:M
    tof = tof_list(ch);
    G_source(:,ch) = [G_source(tof:end, ch) ; zeros(tof-1,1)];
end
tof_list = ones(size(tof_list));


if strcmp(room_memo, 'SALtrunc') == 0
    error('Something got messed up')
end

if rir_database.fs ~= fs
    error('Somethings wrong')
end

%% Optional: Synthetic RIRs

if replace_SAL_with_synth

% TEMP replace SAL with synthetic RIRs
G_source = randn(size(G_source));
room_memo = "synthetic";
rir_desc  = "Synthetic Exponentially Decaying Gaussian";

% Synthetic time delay
synth_delays = 0:(M-1);
synth_delays = synth_delays .* 0;
for ch = 1:M
    delay = synth_delays(ch);
    G_source(:,ch) = [zeros(delay, 1) ; G_source(1:(end - delay),ch)];
end
tof_list = synth_delays + 1;

figure()
for ch = 1:M
    subplot(M,1,ch)
    plot(G_source(:,ch))
    xlabel('Time [samples]')
    title(sprintf("Synthetic RIR Channel %d", ch))
    xlim([0 40])
end

%saveas(gcf, sprintf('%s/Full_T60_Synthetic_RIR.fig', results_dir));


end

%% Optional: Reduce Magnitude of Early reflections
% To make reverberant effect stronger

if enable_reduction_of_ERs


% Manual processing of RIR to make reverberation stronger:
% Make direct sound weaker
for ch = 1:M
    g_SAL_pre    = G_source(:,ch);
    G_source(1:70,ch) = G_source(1:70,ch) ./ 4;

    % EDC_pre  = EDC(g_SAL_pre);
    % EDC_post = EDC(G_source(:,ch));
    
    figure()
    subplot(2,1,1)
    plot((0:(length(g_SAL_pre)-1)) .* (1/fs), g_SAL_pre)
    title(sprintf('Channel %d SAL RIR Before Making early reflections weaker', ch))
    subplot(2,1,2)
    plot((0:(length(G_source(:,ch))-1)) .* (1/fs), G_source(:,ch))
    title(sprintf('Channel %d SAL RIR After Making early reflections weaker', ch))
    
    % figure
    % subplot(2,1,1)
    % plot((0:(length(EDC_pre)-1)) .* (1/fs), 10*log10(EDC_pre));
    % grid on;
    % xlabel('Time [sec]')
    % ylabel('dB')
    % title('EDC for SAL RIR Before Making early reflections weaker')
    % ylim([-65 6])
    % subplot(2,1,2)
    % plot((0:(length(EDC_post)-1)) .* (1/fs), 10*log10(EDC_post));
    % grid on;
    % xlabel('Time [sec]')
    % ylabel('dB')
    % title('EDC for SAL RIR After Making early reflections weaker')
    % ylim([-65 6])

end


end

%% Spatial Noise Data

% Spatial noise
noise_num = 1; % Ventillation office room
S_noise    = rir_database.noise_data.noise_list{noise_num}; 
noise_memo = rir_database.noise_data.noise_memos{noise_num};

% Uncorrelated White Noise
%S_noise = randn(length(s_source), M);
%noise_memo="whiteNoise";

% Loop if needed to synthetically match length to source
num_loops_n = ceil(length(s_source) / length(S_noise(:,1)));
S_noise = repmat(S_noise, num_loops_n, 1);
S_noise = S_noise(1:length(s_source), :);

%% Create results folder

% Results can be saved in the Results Folder
results_dir = sprintf("FINAL_RESULTS/EVAL_T60Sweep_SyntheticRIR_%.0fdB_SNR_%s_%s", SNR_dB, noise_memo, datetime("today"));
run_num = 1;
while exist(results_dir, 'Dir')
    results_dir = sprintf("FINAL_RESULTS/EVAL_T60Sweep_SyntheticRIR_%.0fdB_SNR_%s_%s_%d", SNR_dB, noise_memo, datetime("today"), run_num);
    run_num = run_num + 1;
end
mkdir(results_dir)

%% Design Hearing Aid Gains

fs_h_ha            = 24000;
nfir_24k           = 140;
h_ha_nocomp_24k    = zeros(nfir_24k, 1);
h_ha_nocomp_24k(1) = 1;

HL_freqs_design   = [ 0  (HL_freqs ./ (fs_h_ha/2))  1          ];
HL_HI_design      = [ 0  HL_HI                      HL_HI(end) ];
HL_HI_design      = 10 .^ (HL_HI_design ./ 20);
h_ha_mirror_24k   = fir2(nfir_24k, HL_freqs_design, HL_HI_design);

[h_ha_nalr_24k,~] = eb_NALR(HL_HI, nfir_24k, fs_h_ha); %Design the NAL-R filter

[H_ha_nocomp_24k, freqs_freqz] = freqz(h_ha_nocomp_24k, 1, nfir_24k, fs_h_ha);
[H_ha_mirror_24k,           ~] = freqz(h_ha_mirror_24k, 1, nfir_24k, fs_h_ha);
[H_ha_nalr_24k,             ~] = freqz(h_ha_nalr_24k, 1, nfir_24k, fs_h_ha);

% Design filter halfway between NAL-R and audiogram mirror
H_ha_hybrid_24k_design_dB   = (20*log10(abs(H_ha_nalr_24k)) + 20*log10(abs(H_ha_mirror_24k))) ./ 2;
H_ha_hybrid_24k_design      = 10 .^ (H_ha_hybrid_24k_design_dB/20);
freqs_hybrid_24k_design     = [freqs_freqz'             (fs_h_ha/2) ];
H_ha_hybrid_24k_design      = [H_ha_hybrid_24k_design'  H_ha_hybrid_24k_design(end)];
freqs_hybrid_24k_design     = freqs_hybrid_24k_design ./ (fs_h_ha/2);
freqs_hybrid_24k_design_idx = zeros(length(HL_freqs_design),1);
delta_f = (1 / nfir_24k) / 2;
for f_idx = 1:length(HL_freqs_design)
    %find(abs(freqs_hybrid_24k_design - HL_freqs_design(f_idx)) <= delta_f)
    freqs_hybrid_24k_design_idx(f_idx) = ...
        find(abs(freqs_hybrid_24k_design - HL_freqs_design(f_idx)) <= delta_f);
end
H_ha_hybrid_24k_design      = H_ha_hybrid_24k_design(freqs_hybrid_24k_design_idx);
h_ha_hybrid_24k             = fir2(nfir_24k, HL_freqs_design, H_ha_hybrid_24k_design);
[H_ha_hybrid_24k, ~]        = freqz(h_ha_hybrid_24k, 1, nfir_24k, fs_h_ha);


figure()
plot(HL_freqs, HL_HI, ':k');
hold on;
plot(freqs_freqz, 20*log10(abs(H_ha_nocomp_24k)));
plot(freqs_freqz, 20*log10(abs(H_ha_mirror_24k)));
plot(freqs_freqz, 20*log10(abs(H_ha_hybrid_24k)));
plot(freqs_freqz, 20*log10(abs(H_ha_nalr_24k)));
xlabel('Frequency [Hz]')
ylabel('dB')
legend('True Audiogram Mirror', 'No Compensation', 'Audiogram Mirror Design', 'Hybrid', 'NAL-R Filter')
title('Hearing Aid Gain Designs')

saveas(gcf, sprintf('%s/HearingAidGains.fig', results_dir));

%% Compute SI Plot scale
% Use source * early reflections as the "ideal" test signal and use this to scale NSIM and STMI plots

ER_length   = round((50 / 1000) * fs); % first 50 msec are selected as "early reflections"
ER_length   = min(ER_length, length(G_source(:,1)));
G_source_ER = G_source(1:ER_length, :);

source_data.signal   = s_source;
source_data.rir_data = G_source_ER;
source_data.memo     = source_memo;
source_data.rir_memo = "SAL_ERs";
source_data.rir_desc = "SAL Early Reflections";
source_data.stimdb   = stimdb;

noise_data.enable    = false;
noise_data.signals   = S_noise;
noise_data.SNR_dB    = Inf;
noise_data.memo      = "no noise";

interf_data.enable   = 0;
interf_data.signal   = [0];
interf_data.rir_data = [0];
interf_data.SIR_dB   = Inf;
interf_data.memo     = "None";

[Y_ideal, refstim_ideal] = compute_mic_signals(source_data, noise_data, interf_data, fs);

teststim_ideal  = Y_ideal(:,1);

[metrics_ideal_NH, ~] = run_all_metrics(refstim_ideal, teststim_ideal, fs, HL_NH, stimdb, h_ha_nocomp_24k);

scale_STMI    = 1 / metrics_ideal_NH.STMI;
scale_NSIM_FT = 1 / metrics_ideal_NH.NSIM_FT;
scale_NSIM_MR = 1 / metrics_ideal_NH.NSIM_MR;
scale_HASPI   = 1; % No scaling (already normalized to 0-1)
scale_STOI    = 1; % No scaling (already normalized to 0-1)

plot_scaling_data.metrics_ideal_NH = metrics_ideal_NH;
plot_scaling_data.scale_STMI    = scale_STMI;
plot_scaling_data.scale_NSIM_FT = scale_NSIM_FT;
plot_scaling_data.scale_NSIM_MR = scale_NSIM_MR;
plot_scaling_data.scale_HASPI   = scale_HASPI;
plot_scaling_data.scale_STOI    = scale_STOI;

save(sprintf("%s/plot_scaling_data.mat", results_dir), 'plot_scaling_data');

%% T60 Loop

T60_list = [0 100 250:250:2000] ./ 1000;

STMI_list_NH_unproc    = zeros(length(T60_list),1); 
NSIM_FT_list_NH_unproc = zeros(length(T60_list),1); 
NSIM_MR_list_NH_unproc = zeros(length(T60_list),1); 
HASPI_list_NH_unproc   = zeros(length(T60_list),1); 
STOI_list_NH_unproc    = zeros(length(T60_list),1); 
HASQI_list_NH_unproc   = zeros(length(T60_list),1); 
VISQOL_list_NH_unproc  = zeros(length(T60_list),1); 

STMI_list_HI_unproc    = zeros(length(T60_list),1); 
NSIM_FT_list_HI_unproc = zeros(length(T60_list),1); 
NSIM_MR_list_HI_unproc = zeros(length(T60_list),1); 
HASPI_list_HI_unproc   = zeros(length(T60_list),1); 
STOI_list_HI_unproc    = zeros(length(T60_list),1); 
HASQI_list_HI_unproc   = zeros(length(T60_list),1); 
VISQOL_list_HI_unproc  = zeros(length(T60_list),1); 

STMI_list_NH_proc    = zeros(length(T60_list),1); 
NSIM_FT_list_NH_proc = zeros(length(T60_list),1); 
NSIM_MR_list_NH_proc = zeros(length(T60_list),1); 
HASPI_list_NH_proc   = zeros(length(T60_list),1); 
STOI_list_NH_proc    = zeros(length(T60_list),1); 
HASQI_list_NH_proc   = zeros(length(T60_list),1); 
VISQOL_list_NH_proc  = zeros(length(T60_list),1); 

STMI_list_HI_proc    = zeros(length(T60_list),1); 
NSIM_FT_list_HI_proc = zeros(length(T60_list),1); 
NSIM_MR_list_HI_proc = zeros(length(T60_list),1); 
HASPI_list_HI_proc   = zeros(length(T60_list),1); 
STOI_list_HI_proc    = zeros(length(T60_list),1); 
HASQI_list_HI_proc   = zeros(length(T60_list),1); 
VISQOL_list_HI_proc  = zeros(length(T60_list),1); 

C50_list_unproc  = zeros(length(T60_list),1); 
C50_list_proc    = zeros(length(T60_list),1); 

% Compute mean EDC for SAL room
EDC_list_SAL = zeros(size(G_source));
EDC_SAL      = zeros(length(G_source(:,1)), 1);
for ch = 1:M
    EDC_list_SAL(:, ch) = EDC(G_source(:,ch));
    EDC_SAL = EDC_SAL + EDC_list_SAL(:, ch);
end
EDC_SAL_dB = 10*log10(EDC_SAL / M);

for T60_num = 1:length(T60_list)
    T60 = T60_list(T60_num);
    N60 = T60 * fs;
    fprintf("\n====================================\n")
    fprintf("T60 = %.0f msec\n", T60*1000)
    fprintf("====================================\n")
    
    
    % Channel        
    if T60 == 0
        G_source_trunc = [ones(1,M) ; zeros(1,M)];
    else

        % Compute Additional attenuation needed to get desired T60
        SAL_atten_at_T60 = -1*EDC_SAL_dB(T60*fs);
        added_atten_at_T60 = 60 - SAL_atten_at_T60;
        added_atten_at_T60 = max(added_atten_at_T60, 0);
    
        % Exponential decay curve
        N60       = T60 * fs;
        %L_channel = min(round(N60*1.25), length(G_source(:,1)));
        L_channel = min(round(N60*2), length(G_source(:,1)));
        tau       = N60 / log(10^(added_atten_at_T60 / 20)); % -ln(10^(-60/20)) = ln(10^3)
        
        % Apply synthetic exponential decay to 2.1sec MYRiAD real RIR
        G_source_trunc = G_source(1:L_channel, :);
        %figure()
        for ch = 1:M
            tof = tof_list(ch);
            exp_decay = [ones(tof-1,1) ; exp(-1 .* (0:(L_channel-tof))' ./ tau)];
            G_source_trunc(:,ch)  = G_source_trunc(:,ch)  .* exp_decay;   

            g_SAL_ch      = G_source(:,ch);
            g_SALtrunc_ch = G_source_trunc(:,ch);
            
            edc_SAL_ch      = EDC(g_SAL_ch);
            edc_exp          = EDC(exp_decay);
            edc_SALtrunc_ch4 = EDC(g_SALtrunc_ch);
            
            % subplot(M, 1, ch)
            % plot((0:(length(edc_SALtrunc_ch4)-1)) .* (1/fs), 10*log10(edc_SALtrunc_ch4));
            % grid on;
            % xlabel('Time [sec]')
            % ylabel('dB')
            % plt_title = sprintf("EDC for Truncated SAL RIR (T60 = %.0f msec)", T60*1000);
            % title(plt_title)
            % ylim([-80 6])
            % xlim([0 3])

        end

        figure()
        for ch = 1:M
            subplot(M,1,ch)
            plot(G_source_trunc(:,ch))
            xlabel('Time [samples]')
            title(sprintf("RIR Channel %d", ch))
            xlim([1 N60])
        end
    end


    if enable_debug_plots == 1


        g_SAL_ch      = G_source(:,4);
        g_SALtrunc_ch = G_source_trunc(:,4);
    
        edc_SAL_ch      = EDC(g_SAL_ch);
        edc_exp          = EDC(exp_decay);
        edc_SALtrunc_ch4 = EDC(g_SALtrunc_ch);
        
        figure()
        subplot(3,1,1)
        plot((0:(length(g_SAL_ch)-1)) .* (1/fs), g_SAL_ch);
        title('SAL RIR (T60 = 2100 msec)')
        xlabel('Time [sec]')
        xlim([0 3])
        subplot(3,1,2)
        plot((0:(length(exp_decay)-1)) .* (1/fs), exp_decay);
        title('Applied Exponential Decay Function')
        xlabel('Time [sec]')
        xlim([0 3])
        subplot(3,1,3)
        plot((0:(length(g_SALtrunc_ch)-1)) .* (1/fs), g_SALtrunc_ch);
        plt_title = sprintf("Truncated SAL RIR (T60 = %.0f msec)", T60 * 1000);
        title(plt_title)
        xlabel('Time [sec]')
        xlim([0 3])
    
        figure()
        subplot(3,1,1)
        plot((0:(length(edc_SAL_ch)-1)) .* (1/fs), 10*log10(edc_SAL_ch));
        grid on;
        xlabel('Time [sec]')
        ylabel('dB')
        title('EDC for SAL RIR (T60 = 2100 msec)')
        ylim([-80 6])
        xlim([0 3])
        subplot(3,1,2)
        plot((0:(length(edc_exp)-1)) .* (1/fs), 10*log10(edc_exp));
        grid on;
        xlabel('Time [sec]')
        ylabel('dB')
        title('EDC for Applied Exponential Decay Function')
        ylim([-80 6])
        xlim([0 3])
        subplot(3,1,3)
        plot((0:(length(edc_SALtrunc_ch4)-1)) .* (1/fs), 10*log10(edc_SALtrunc_ch4));
        grid on;
        xlabel('Time [sec]')
        ylabel('dB')
        plt_title = sprintf("EDC for Truncated SAL RIR (T60 = %.0f msec)", T60*1000);
        title(plt_title)
        ylim([-80 6])
        xlim([0 3])

    end
    
    source_data.signal   = s_source;
    source_data.rir_data = G_source_trunc;
    source_data.memo     = source_memo;
    source_data.rir_memo = room_memo;
    source_data.rir_desc = rir_desc;
    source_data.stimdb   = stimdb;
    
    noise_data.enable    = enable_noise;
    noise_data.signals   = S_noise;
    noise_data.SNR_dB    = SNR_dB;
    noise_data.memo      = noise_memo;
    
    interf_data.enable   = 0;
    interf_data.signal   = [0];
    interf_data.rir_data = [0];
    interf_data.SIR_dB   = Inf;
    interf_data.memo     = "None";
    
    [metrics] = run_eval(source_data, ...
                         noise_data,  ...
                         interf_data, ...
                         fs, ...
                         N60, ...
                         s1_enable, ...
                         s1_on_clean_speech, ...
                         HL_NH, ...
                         HL_HI, ...
                         h_ha_nocomp_24k, ...
                         h_ha_nalr_24k, ...
                         results_dir);

    % Save to buffers
    STMI_list_NH_unproc(T60_num)    = metrics.metrics_unproc_NH.STMI;
    NSIM_FT_list_NH_unproc(T60_num) = metrics.metrics_unproc_NH.NSIM_FT;
    NSIM_MR_list_NH_unproc(T60_num) = metrics.metrics_unproc_NH.NSIM_MR;
    HASPI_list_NH_unproc(T60_num)   = metrics.metrics_unproc_NH.HASPI;
    STOI_list_NH_unproc(T60_num)    = metrics.metrics_unproc_NH.STOI;
    HASQI_list_NH_unproc(T60_num)   = metrics.metrics_unproc_NH.HASQI;
    VISQOL_list_NH_unproc(T60_num)  = metrics.metrics_unproc_NH.VISQOL;

    STMI_list_HI_unproc(T60_num)    = metrics.metrics_unproc_HI.STMI;
    NSIM_FT_list_HI_unproc(T60_num) = metrics.metrics_unproc_HI.NSIM_FT;
    NSIM_MR_list_HI_unproc(T60_num) = metrics.metrics_unproc_HI.NSIM_MR;
    HASPI_list_HI_unproc(T60_num)   = metrics.metrics_unproc_HI.HASPI;
    STOI_list_HI_unproc(T60_num)    = metrics.metrics_unproc_HI.STOI;
    HASQI_list_HI_unproc(T60_num)   = metrics.metrics_unproc_HI.HASQI;
    VISQOL_list_HI_unproc(T60_num)  = metrics.metrics_unproc_HI.VISQOL;

    STMI_list_NH_proc(T60_num)    = metrics.metrics_proc_NH.STMI;
    NSIM_FT_list_NH_proc(T60_num) = metrics.metrics_proc_NH.NSIM_FT;
    NSIM_MR_list_NH_proc(T60_num) = metrics.metrics_proc_NH.NSIM_MR;
    HASPI_list_NH_proc(T60_num)   = metrics.metrics_proc_NH.HASPI;
    STOI_list_NH_proc(T60_num)    = metrics.metrics_proc_NH.STOI;
    HASQI_list_NH_proc(T60_num)   = metrics.metrics_proc_NH.HASQI;
    VISQOL_list_NH_proc(T60_num)  = metrics.metrics_proc_NH.VISQOL;

    STMI_list_HI_proc(T60_num)    = metrics.metrics_proc_HI.STMI;
    NSIM_FT_list_HI_proc(T60_num) = metrics.metrics_proc_HI.NSIM_FT;
    NSIM_MR_list_HI_proc(T60_num) = metrics.metrics_proc_HI.NSIM_MR;
    HASPI_list_HI_proc(T60_num)   = metrics.metrics_proc_HI.HASPI;
    STOI_list_HI_proc(T60_num)    = metrics.metrics_proc_HI.STOI;
    HASQI_list_HI_proc(T60_num)   = metrics.metrics_proc_HI.HASQI;
    VISQOL_list_HI_proc(T60_num)  = metrics.metrics_proc_HI.VISQOL;

    C50_list_unproc(T60_num) = metrics.metrics_physical_unproc.C50;
    C50_list_proc(T60_num)   = metrics.metrics_physical_proc.C50;
     
end

%% PLOTS

test_memo_NH   = sprintf("(SAL RIR truncated, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f], No HA Gain)", SNR_dB, stimdb, HL_NH(1), HL_NH(2), HL_NH(3), HL_NH(4), HL_NH(5), HL_NH(6));
test_memo_HI   = sprintf("(SAL RIR truncated, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f], NAL-R Gain)", SNR_dB, stimdb, HL_HI(1), HL_HI(2), HL_HI(3), HL_HI(4), HL_HI(5), HL_HI(6));
test_memo_phys = sprintf("(SAL RIR truncated, SNR = %.0f dB, Stimulus = %.0f dBSPL)", SNR_dB, stimdb);

plot_colours = get(gca,'colororder');

% SPEECH INTELLIGIBILITY vs T60

figure()
set(gcf,'Position',fig_size)
subplot(1,2,1)
plot(T60_list, STMI_list_NH_unproc,    'Color', plot_colours(1,:), 'LineStyle', '--', 'Marker', 'o')
hold on;
plot(T60_list, STMI_list_NH_proc,      'Color', plot_colours(1,:), 'LineStyle', '-',  'Marker', 'o')
plot(T60_list, NSIM_FT_list_NH_unproc, 'Color', plot_colours(2,:), 'LineStyle', '--', 'Marker', '^')
plot(T60_list, NSIM_FT_list_NH_proc,   'Color', plot_colours(2,:), 'LineStyle', '-',  'Marker', '^')
plot(T60_list, NSIM_MR_list_NH_unproc, 'Color', plot_colours(3,:), 'LineStyle', '--', 'Marker', 'diamond')
plot(T60_list, NSIM_MR_list_NH_proc,   'Color', plot_colours(3,:), 'LineStyle', '-',  'Marker', 'diamond')
plot(T60_list, HASPI_list_NH_unproc,   'Color', plot_colours(4,:), 'LineStyle', '--', 'Marker', 'square')
plot(T60_list, HASPI_list_NH_proc,     'Color', plot_colours(4,:), 'LineStyle', '-',  'Marker', 'square')
plot(T60_list, STOI_list_NH_unproc,    'Color', plot_colours(5,:), 'LineStyle', '--', 'Marker', '+')
plot(T60_list, STOI_list_NH_proc,      'Color', plot_colours(5,:), 'LineStyle', '-',  'Marker', '+')
ylim([0 Inf])
xlabel('T60 (sec)')
legend("STMI (unprocessed)",    "STMI (processed)", ...
       "NSIM FT (unprocessed)", "NSIM FT (processed)", ...
       "NSIM MR (unprocessed)", "NSIM MR (processed)", ...
       "HASPI (unprocessed)",   "HASPI (processed)", ...
       "STOI (unprocessed)",    "STOI (processed)")
title("Speech Intelligibility Predictors", test_memo_NH)

subplot(1,2,2)
plot(T60_list, STMI_list_HI_unproc,    'Color', plot_colours(1,:), 'LineStyle', '--', 'Marker', 'o')
hold on;
plot(T60_list, STMI_list_HI_proc,      'Color', plot_colours(1,:), 'LineStyle', '-',  'Marker', 'o')
plot(T60_list, NSIM_FT_list_HI_unproc, 'Color', plot_colours(2,:), 'LineStyle', '--', 'Marker', '^')
plot(T60_list, NSIM_FT_list_HI_proc,   'Color', plot_colours(2,:), 'LineStyle', '-',  'Marker', '^')
plot(T60_list, NSIM_MR_list_HI_unproc, 'Color', plot_colours(3,:), 'LineStyle', '--', 'Marker', 'diamond')
plot(T60_list, NSIM_MR_list_HI_proc,   'Color', plot_colours(3,:), 'LineStyle', '-',  'Marker', 'diamond')
plot(T60_list, HASPI_list_HI_unproc,   'Color', plot_colours(4,:), 'LineStyle', '--', 'Marker', 'square')
plot(T60_list, HASPI_list_HI_proc,     'Color', plot_colours(4,:), 'LineStyle', '-',  'Marker', 'square')
ylim([0 Inf])
xlabel('T60 (sec)')
legend("STMI (unprocessed)",    "STMI (processed)", ...
       "NSIM FT (unprocessed)", "NSIM FT (processed)", ...
       "NSIM MR (unprocessed)", "NSIM MR (processed)", ...
       "HASPI (unprocessed)",   "HASPI (processed)")
ylim([0 Inf])
xlabel('T60 (sec)')
title("Speech Intelligibility Predictors", test_memo_HI)

saveas(gcf, sprintf('%s/SI_v_T60.fig', results_dir));

figure()
set(gcf,'Position',fig_size)
subplot(1,2,1)
plot(T60_list, STMI_list_NH_unproc .* scale_STMI,    'Color', plot_colours(1,:), 'LineStyle', '--', 'Marker', 'o')
hold on;
plot(T60_list, STMI_list_NH_proc .* scale_STMI,      'Color', plot_colours(1,:), 'LineStyle', '-',  'Marker', 'o')
plot(T60_list, NSIM_FT_list_NH_unproc .* scale_NSIM_FT, 'Color', plot_colours(2,:), 'LineStyle', '--', 'Marker', '^')
plot(T60_list, NSIM_FT_list_NH_proc .* scale_NSIM_FT,   'Color', plot_colours(2,:), 'LineStyle', '-',  'Marker', '^')
plot(T60_list, NSIM_MR_list_NH_unproc .* scale_NSIM_MR, 'Color', plot_colours(3,:), 'LineStyle', '--', 'Marker', 'diamond')
plot(T60_list, NSIM_MR_list_NH_proc .* scale_NSIM_MR,   'Color', plot_colours(3,:), 'LineStyle', '-',  'Marker', 'diamond')
plot(T60_list, HASPI_list_NH_unproc .* scale_HASPI,   'Color', plot_colours(4,:), 'LineStyle', '--', 'Marker', 'square')
plot(T60_list, HASPI_list_NH_proc .* scale_HASPI,     'Color', plot_colours(4,:), 'LineStyle', '-',  'Marker', 'square')
plot(T60_list, STOI_list_NH_unproc .* scale_STOI,    'Color', plot_colours(5,:), 'LineStyle', '--', 'Marker', '+')
plot(T60_list, STOI_list_NH_proc .* scale_STOI,      'Color', plot_colours(5,:), 'LineStyle', '-',  'Marker', '+')
ylim([0 Inf])
xlabel('T60 (sec)')
legend("STMI (unprocessed)",    "STMI (processed)", ...
       "NSIM FT (unprocessed)", "NSIM FT (processed)", ...
       "NSIM MR (unprocessed)", "NSIM MR (processed)", ...
       "HASPI (unprocessed)",   "HASPI (processed)", ...
       "STOI (unprocessed)",    "STOI (processed)")
title("Scaled Speech Intelligibility Predictors", test_memo_NH)

subplot(1,2,2)
plot(T60_list, STMI_list_HI_unproc .* scale_STMI,    'Color', plot_colours(1,:), 'LineStyle', '--', 'Marker', 'o')
hold on;
plot(T60_list, STMI_list_HI_proc .* scale_STMI,      'Color', plot_colours(1,:), 'LineStyle', '-',  'Marker', 'o')
plot(T60_list, NSIM_FT_list_HI_unproc .* scale_NSIM_FT, 'Color', plot_colours(2,:), 'LineStyle', '--', 'Marker', '^')
plot(T60_list, NSIM_FT_list_HI_proc .* scale_NSIM_FT,   'Color', plot_colours(2,:), 'LineStyle', '-',  'Marker', '^')
plot(T60_list, NSIM_MR_list_HI_unproc .* scale_NSIM_MR, 'Color', plot_colours(3,:), 'LineStyle', '--', 'Marker', 'diamond')
plot(T60_list, NSIM_MR_list_HI_proc .* scale_NSIM_MR,   'Color', plot_colours(3,:), 'LineStyle', '-',  'Marker', 'diamond')
plot(T60_list, HASPI_list_HI_unproc .* scale_HASPI,   'Color', plot_colours(4,:), 'LineStyle', '--', 'Marker', 'square')
plot(T60_list, HASPI_list_HI_proc .* scale_HASPI,     'Color', plot_colours(4,:), 'LineStyle', '-',  'Marker', 'square')
ylim([0 Inf])
xlabel('T60 (sec)')
legend("STMI (unprocessed)",    "STMI (processed)", ...
       "NSIM FT (unprocessed)", "NSIM FT (processed)", ...
       "NSIM MR (unprocessed)", "NSIM MR (processed)", ...
       "HASPI (unprocessed)",   "HASPI (processed)")
ylim([0 Inf])
xlabel('T60 (sec)')
title("Scaled Speech Intelligibility Predictors", test_memo_HI)

saveas(gcf, sprintf('%s/SI_v_T60_scaled.fig', results_dir));

% SPEECH QUALITY vs T60

max_SQ = max([max(HASQI_list_NH_unproc) max(HASQI_list_NH_proc) max(VISQOL_list_NH_unproc) max(VISQOL_list_NH_proc) max(HASQI_list_HI_unproc) max(HASQI_list_HI_proc)]);
ylim_SQ = max_SQ*1.1;

figure()
set(gcf,'Position',fig_size)
subplot(1,2,1)
plot(T60_list, HASQI_list_NH_unproc,  'Color', plot_colours(4,:), 'LineStyle', '--', 'Marker', 'square')
hold on;
plot(T60_list, HASQI_list_NH_proc,    'Color', plot_colours(4,:), 'LineStyle', '-',  'Marker', 'square')
plot(T60_list, VISQOL_list_NH_unproc, 'Color', plot_colours(1,:), 'LineStyle', '--', 'Marker', 'diamond')
plot(T60_list, VISQOL_list_NH_proc,   'Color', plot_colours(1,:), 'LineStyle', '-',  'Marker', 'diamond')
ylim([0 ylim_SQ])

xlabel('T60 (sec)')
legend("HASQI (unprocessed)",  "HASQI (processed)", ...
       "VISQOL (unprocessed)", "VISQOL (processed)")
title("Speech Quality Predictors", test_memo_NH)

subplot(1,2,2)
plot(T60_list, HASQI_list_HI_unproc,  'Color', plot_colours(4,:), 'LineStyle', '--', 'Marker', 'square')
hold on;
plot(T60_list, HASQI_list_HI_proc,    'Color', plot_colours(4,:), 'LineStyle', '-',  'Marker', 'square')
xlabel('T60 (sec)')
legend("HASQI (unprocessed)", "HASQI (processed)")
title("Speech Quality Predictors", test_memo_HI)
ylim([0 ylim_SQ])

saveas(gcf, sprintf('%s/SQ_v_T60.fig', results_dir));

% CLARITY vs T60

figure()
%set(gcf,'Position',fig_size)
plot(T60_list, 10*log10(C50_list_unproc), "--square")
hold on;
plot(T60_list, 10*log10(C50_list_proc), "-square")
xlabel('T60 (sec)')
ylabel('dB')
legend("C50 (unprocessed)", "C50 (processed)")
title("Clarity (C50)", test_memo_phys)

saveas(gcf, sprintf('%s/C50_v_T60.fig', results_dir));

save(sprintf("%s/STMI_list_NH_unproc.mat", results_dir), 'STMI_list_NH_unproc');
save(sprintf("%s/NSIM_FT_list_NH_unproc.mat", results_dir), 'NSIM_FT_list_NH_unproc');
save(sprintf("%s/NSIM_MR_list_NH_unproc.mat", results_dir), 'NSIM_MR_list_NH_unproc');
save(sprintf("%s/HASPI_list_NH_unproc.mat", results_dir), 'HASPI_list_NH_unproc');
save(sprintf("%s/STOI_list_NH_unproc.mat", results_dir), 'STOI_list_NH_unproc');
save(sprintf("%s/HASQI_list_NH_unproc.mat", results_dir), 'HASQI_list_NH_unproc');
save(sprintf("%s/VISQOL_list_NH_unproc.mat", results_dir), 'VISQOL_list_NH_unproc');

save(sprintf("%s/STMI_list_NH_proc.mat", results_dir), 'STMI_list_NH_proc');
save(sprintf("%s/NSIM_FT_list_NH_proc.mat", results_dir), 'NSIM_FT_list_NH_proc');
save(sprintf("%s/NSIM_MR_list_NH_proc.mat", results_dir), 'NSIM_MR_list_NH_proc');
save(sprintf("%s/HASPI_list_NH_proc.mat", results_dir), 'HASPI_list_NH_proc');
save(sprintf("%s/STOI_list_NH_proc.mat", results_dir), 'STOI_list_NH_proc');
save(sprintf("%s/HASQI_list_NH_proc.mat", results_dir), 'HASQI_list_NH_proc');
save(sprintf("%s/VISQOL_list_NH_proc.mat", results_dir), 'VISQOL_list_NH_proc');

save(sprintf("%s/STMI_list_HI_unproc.mat", results_dir), 'STMI_list_HI_unproc');
save(sprintf("%s/NSIM_FT_list_HI_unproc.mat", results_dir), 'NSIM_FT_list_HI_unproc');
save(sprintf("%s/NSIM_MR_list_HI_unproc.mat", results_dir), 'NSIM_MR_list_HI_unproc');
save(sprintf("%s/HASPI_list_HI_unproc.mat", results_dir), 'HASPI_list_HI_unproc');
save(sprintf("%s/HASQI_list_HI_unproc.mat", results_dir), 'HASQI_list_HI_unproc');

save(sprintf("%s/STMI_list_HI_proc.mat", results_dir), 'STMI_list_HI_proc');
save(sprintf("%s/NSIM_FT_list_HI_proc.mat", results_dir), 'NSIM_FT_list_HI_proc');
save(sprintf("%s/NSIM_MR_list_HI_proc.mat", results_dir), 'NSIM_MR_list_HI_proc');
save(sprintf("%s/HASPI_list_HI_proc.mat", results_dir), 'HASPI_list_HI_proc');
save(sprintf("%s/HASQI_list_HI_proc.mat", results_dir), 'HASQI_list_HI_proc');

save(sprintf("%s/C50_list_unproc.mat", results_dir), 'C50_list_unproc');
save(sprintf("%s/C50_list_proc.mat", results_dir), 'C50_list_proc');
