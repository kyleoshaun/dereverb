restoredefaultpath

close all
clear
clc

% myCluster = parcluster('Processes');
% delete(myCluster.Jobs)
% delete(gcp('nocreate'))
% pause(5)
% parpool

%% Paths

addpath('../dereverb/LPC/../../samples')
addpath('../dereverb/LPC/../../RIR_Databases/AIR_1_4_BinauralDatabase/')
addpath('../dereverb/LPC/../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/')
addpath('../dereverb/LPC/../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/')
addpath('../dereverb/LPC/../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/courtyard/')

% Integration of EC with BEZ model
addpath("./BSIM2020_EC_BEZ2018a_model/")

% HASPI/HASQI
addpath("./HASPIv2_HASQIv2_HAAQIv1_Common/")

% BSIM
addpath('BSIM_2020/')
addpath('BSIM_2020')
addpath('BSIM_2020/SRMRToolbox-master')
addpath('BSIM_2020/SRMRToolbox-master/libs/')
addpath('BSIM_2020/SRMRToolbox-master/libs/vad/')
addpath('BSIM_2020/SRMRToolbox-master/libs/PreProc/')
addpath('BSIM_2020/SRMRToolbox-master/libs/gammatonegram/')
addpath('BSIM_2020/SRMRToolbox-master/libs/auditory/')

addpath(genpath('BSIM_2020/Anechoic/'))
addpath('BSIM_2020/functions/')
addpath('BSIM_2020/functions/gammatonefilter/')
addpath('BSIM_2020/functions/ltfat/')

% Model (BEZ2018a)
addpath('BEZ2018a_model/')

% Other
addpath('NSIM_Papers_ExampleCode 2/NSIMforBEZ2018a_model/')
addpath('STMI_Papers_ExampleCode 1/STMIforBEZ2018a_model')

addpath("../dereverb/utilities/matlab/")

%% Shuffle and save seed for reproducibility
rng('shuffle')
seed = abs(round(randn(1,1) * 10000));
rng(seed);
rng(724); % Reproduce specific seed


%% Init

fs = 16000;

SNR_dB = 300;

HL_freqs = [250 500 1000 2000 4000 6000]; % 6 audiometric frequencies used by HASPI
HL_NH = [0   0   0  0  0  0]; % Normal Hearing
HL_HI = [35  35  40 50 60 65]; % High-freq hearing loss (IEC 60118-15 Moderate HL, Moderately Sloping Group)

[sref_precalib, Fs_stim] = audioread("SA1.WAV");

if Fs_stim ~= fs
    sref_precalib = resample(sref_precalib, fs, Fs_stim);
end

stimdb = 65; % speech level in dB SPL 

fig_size = [10 10 1300 600];


%% Load real RIR data

rir_database = load('../RIR_Databases/rir_databases_4micCombined_rirOnly.mat');
rir_database = rir_database.rir_databases_4micCombined_rir_only;

if rir_database.fs ~= fs
    error("RIR Database sample rate mismatch");
end

room_num       = 4; % SAL (T60 = 2.1 sec)
G_source       = rir_database.rir_data.rir_list_0deg{room_num};
tof_list       = rir_database.rir_data.tof_est{room_num};
room_memo      = rir_database.rir_data.rir_memos{room_num};

if strcmp(room_memo, 'SAL') == 0
    error('Something got messed up')
end

if rir_database.fs ~= fs
    error('Somethings wrong')
end

g_left_SAL  = G_source(:,1);
g_right_SAL = G_source(:,2);

% Manual processing of RIR to make reverberation stronger:
% Make direct sound weaker
g_left_SAL_pre   = g_left_SAL;
g_left_SAL(1:70) = g_left_SAL(1:70) ./ 2;

EDC_pre  = EDC(g_left_SAL_pre);
EDc_post = EDC(g_left_SAL);

figure()
subplot(2,1,1)
plot((0:(length(g_left_SAL_pre)-1)) .* (1/fs), g_left_SAL_pre)
title('SAL RIR Before Making early reflections weaker')
subplot(2,1,2)
plot((0:(length(g_left_SAL)-1)) .* (1/fs), g_left_SAL)
title('SAL RIR After Making early reflections weaker')

figure
subplot(2,1,1)
plot((0:(length(EDC_pre)-1)) .* (1/fs), 10*log10(EDC_pre));
grid on;
xlabel('Time [sec]')
ylabel('dB')
title('EDC for SAL RIR Before Making early reflections weaker')
ylim([-65 6])
subplot(2,1,2)
plot((0:(length(EDc_post)-1)) .* (1/fs), 10*log10(EDc_post));
grid on;
xlabel('Time [sec]')
ylabel('dB')
title('EDC for SAL RIR After Making early reflections weaker')
ylim([-65 6])


%% Create results folder

% Results can be saved in the Results Folder
results_dir = sprintf("Results_T60Sweep_TruncatedSALProcessedERs_%.0fdB_SNR_FIXED_dBSPL", SNR_dB);
if ~exist(results_dir,'Dir')
    mkdir(results_dir)
end

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


%% T60 Loop

T60_list = (0:100:2500) ./ 1000;

STMI_list_NH           = zeros(length(T60_list),1); 
NSIM_FT_list_NH        = zeros(length(T60_list),1); 
NSIM_MR_list_NH        = zeros(length(T60_list),1); 
HASPI_list_NH          = zeros(length(T60_list),1); 
STOI_list_NH           = zeros(length(T60_list),1); 
HASQI_list_NH          = zeros(length(T60_list),1); 
VISQOL_list_NH         = zeros(length(T60_list),1); 

STMI_list_HI           = zeros(length(T60_list),1); 
NSIM_FT_list_HI        = zeros(length(T60_list),1); 
NSIM_MR_list_HI        = zeros(length(T60_list),1); 
HASPI_list_HI          = zeros(length(T60_list),1); 
STOI_list_HI           = zeros(length(T60_list),1); 
HASQI_list_HI          = zeros(length(T60_list),1); 
VISQOL_list_HI         = zeros(length(T60_list),1); 

for T60_num = 1:length(T60_list)
    T60 = T60_list(T60_num);
    %fprintf("\n====================================\n")
    fprintf("T60 = %.0f msec\n", T60*1000)
    %fprintf("====================================\n")
    
    
    % Channel
    
    channel_memo = "synthetic (exponentially decaying white noise)";
        
    if T60 == 0
        g_left  = 1;
        g_right = 1;
    else

        % Compute Additional attenuation needed to get desired T60
        SAL_dB_per_sec     = 44.7 / 1.5;
        SAL_atten_at_T60   = SAL_dB_per_sec * T60 + 2;
        added_atten_at_T60 = 60 - SAL_atten_at_T60;
    
        % Exponential decay curve
        N60 = T60 * fs;
        L_channel = length(g_left_SAL);
        tau = N60 / log(10^(added_atten_at_T60 / 20)); % -ln(10^(-60/20)) = ln(10^3)
        exp_decay = exp(-1 .* (0:(L_channel-1))' ./ tau);
        
        % Apply synthetic exponential decay to 2.1sec MYRiAD real RIR
        g_left  = g_left_SAL  .* exp_decay;   
        g_right = g_right_SAL .* exp_decay;   

    end


    % edc_myriad = EDC(g_left_SAL);
    % edc_exp    = EDC(exp_decay);
    % edc_left   = EDC(g_left);
    % figure()
    % subplot(3,1,1)
    % plot((0:(length(g_left_SAL)-1)) .* (1/fs), g_left_SAL);
    % title('SAL RIR (T60 = 2100 msec)')
    % xlabel('Time [sec]')
    % subplot(3,1,2)
    % plot((0:(length(exp_decay)-1)) .* (1/fs), exp_decay);
    % title('Applied Exponential Decay Function')
    % xlabel('Time [sec]')
    % subplot(3,1,3)
    % plot((0:(length(g_left)-1)) .* (1/fs), g_left);
    % plt_title = sprintf("Truncated SAL RIR (T60 = %.0f msec)", T60 * 1000);
    % title(plt_title)
    % xlabel('Time [sec]')
    % 
    % figure()
    % subplot(3,1,1)
    % plot((0:(length(edc_myriad)-1)) .* (1/fs), 10*log10(edc_myriad));
    % grid on;
    % xlabel('Time [sec]')
    % ylabel('dB')
    % title('EDC for SAL RIR (T60 = 2100 msec)')
    % ylim([-80 6])
    % xlim([0 3])
    % subplot(3,1,2)
    % plot((0:(length(edc_myriad)-1)) .* (1/fs), 10*log10(edc_exp));
    % grid on;
    % xlabel('Time [sec]')
    % ylabel('dB')
    % title('EDC for Applied Exponential Decay Function')
    % ylim([-80 6])
    % xlim([0 3])
    % subplot(3,1,3)
    % plot((0:(length(edc_myriad)-1)) .* (1/fs), 10*log10(edc_left));
    % grid on;
    % xlabel('Time [sec]')
    % ylabel('dB')
    % plt_title = sprintf("EDC for Truncated SAL RIR (T60 = %.0f msec)", T60*1000);
    % title(plt_title)
    % ylim([-80 6])
    % xlim([0 3])
    
    % Mic Signals
    teststim  = filter(g_left,  1, sref_precalib);

    % Calibrate dB SPL
    refstim   = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim  = teststim/rms(teststim)*20e-6*10^(stimdb/20);

    % Add noise
    SNR = 10 ^ (SNR_dB / 20);
    s_noise  = randn(length(teststim),   1);
    SNR_i       = rms(teststim) / rms(s_noise);
    teststim  = teststim  + (s_noise ./ (SNR / SNR_i));

    % Run
    [metrics_NH, nalr_struct_NH] = run_all_metrics(refstim, teststim, fs, HL_NH, stimdb, h_ha_nocomp_24k);
    [metrics_HI, nalr_struct_HI] = run_all_metrics(refstim, teststim, fs, HL_HI, stimdb, h_ha_nalr_24k);

    % Save NAL-R Filter info
    audiowrite(sprintf("%s/test_stim_noNALR_T60_%.0fmsec.wav",   results_dir, T60*1000), teststim .* (0.5 / max(abs(teststim))), fs);
    audiowrite(sprintf("%s/test_stim_withNALR_T60_%.0fmsec.wav", results_dir, T60*1000), nalr_struct_HI.teststim_post_nalr .* (0.5 / max(abs(nalr_struct_HI.teststim_post_nalr))), fs);
    save(sprintf("%s/nalr_struct_HI_T60_%.0fmsec.mat", results_dir, T60*1000),"-fromstruct", nalr_struct_HI);

    % Save to buffers
    STMI_list_NH(T60_num)           = metrics_NH.STMI;
    NSIM_FT_list_NH(T60_num)        = metrics_NH.NSIM_FT;
    NSIM_MR_list_NH(T60_num)        = metrics_NH.NSIM_MR;
    HASPI_list_NH(T60_num)          = metrics_NH.HASPI;
    STOI_list_NH(T60_num)           = metrics_NH.STOI;
    HASQI_list_NH(T60_num)          = metrics_NH.HASQI;
    VISQOL_list_NH(T60_num)         = metrics_NH.VISQOL;

    STMI_list_HI(T60_num)           = metrics_HI.STMI;
    NSIM_FT_list_HI(T60_num)        = metrics_HI.NSIM_FT;
    NSIM_MR_list_HI(T60_num)        = metrics_HI.NSIM_MR;
    HASPI_list_HI(T60_num)          = metrics_HI.HASPI;
    STOI_list_HI(T60_num)           = metrics_HI.STOI;
    HASQI_list_HI(T60_num)          = metrics_HI.HASQI;
    VISQOL_list_HI(T60_num)         = metrics_HI.VISQOL;
     
end

test_memo_NH = sprintf("(SAL RIR truncated w/ ERs reduced by 6 dB, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f], No HA Gain)", SNR_dB, stimdb, HL_NH(1), HL_NH(2), HL_NH(3), HL_NH(4), HL_NH(5), HL_NH(6));
test_memo_HI = sprintf("(SAL RIR truncated w/ ERs reduced by 6 dB, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f], NAL-R Gain)", SNR_dB, stimdb, HL_HI(1), HL_HI(2), HL_HI(3), HL_HI(4), HL_HI(5), HL_HI(6));

figure()
set(gcf,'Position',fig_size)
subplot(1,2,1)
plot(T60_list, STMI_list_NH, "-o")
hold on;
plot(T60_list, NSIM_FT_list_NH, "-^")
plot(T60_list, NSIM_MR_list_NH, "-diamond")
plot(T60_list, HASPI_list_NH, "-square")
plot(T60_list, STOI_list_NH, "-+")
ylim([0 Inf])
xlabel('T60 (sec)')
legend("STMI", "NSIM FT", "NSIM MR", "HASPI", "STOI")
title("Speech Intelligibility Predictors", test_memo_NH)

subplot(1,2,2)
plot(T60_list, STMI_list_HI, "-o")
hold on;
plot(T60_list, NSIM_FT_list_HI, "-^")
plot(T60_list, NSIM_MR_list_HI, "-diamond")
plot(T60_list, HASPI_list_HI, "-square")
ylim([0 Inf])
xlabel('T60 (sec)')
legend("STMI", "NSIM FT", "NSIM MR", "HASPI")
title("Speech Intelligibility Predictors", test_memo_HI)

saveas(gcf, sprintf('%s/SI_v_T60.fig', results_dir));

figure()
set(gcf,'Position',fig_size)
subplot(1,2,1)
plot(T60_list, HASQI_list_NH, "-square")
hold on;
plot(T60_list, VISQOL_list_NH, "-diamond")
xlabel('T60 (sec)')
legend("HASQI", "VISQOL")
title("Speech Quality Predictors", test_memo_NH)

subplot(1,2,2)
plot(T60_list, HASQI_list_HI, "-square")
xlabel('T60 (sec)')
legend("HASQI")
title("Speech Quality Predictors", test_memo_HI)

saveas(gcf, sprintf('%s/SQ_v_T60.fig', results_dir));

save(sprintf("%s/STMI_list_NH.mat", results_dir), 'STMI_list_NH');
save(sprintf("%s/NSIM_FT_list_NH.mat", results_dir), 'NSIM_FT_list_NH');
save(sprintf("%s/NSIM_MR_list_NH.mat", results_dir), 'NSIM_MR_list_NH');
save(sprintf("%s/HASPI_list_NH.mat", results_dir), 'HASPI_list_NH');
save(sprintf("%s/STOI_list_NH.mat", results_dir), 'STOI_list_NH');
save(sprintf("%s/HASQI_list_NH.mat", results_dir), 'HASQI_list_NH');
save(sprintf("%s/VISQOL_list_NH.mat", results_dir), 'VISQOL_list_NH');

save(sprintf("%s/STMI_list_HI.mat", results_dir), 'STMI_list_HI');
save(sprintf("%s/NSIM_FT_list_HI.mat", results_dir), 'NSIM_FT_list_HI');
save(sprintf("%s/NSIM_MR_list_HI.mat", results_dir), 'NSIM_MR_list_HI');
save(sprintf("%s/HASPI_list_HI.mat", results_dir), 'HASPI_list_HI');
save(sprintf("%s/HASQI_list_HI.mat", results_dir), 'HASQI_list_HI');
