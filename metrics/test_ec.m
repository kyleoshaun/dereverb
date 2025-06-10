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

addpath('../dereverb/utilities/matlab')

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

%% Shuffle and save seed for reproducibility
rng('shuffle')
seed = abs(round(randn(1,1) * 10000));
rng(seed);
%rng(724); % Reproduce specific seed

%% TEMP THIS SHOULD BE IN EC_FrontEnd.m

Model_init;


%% Init

[sref_precalib, Fs_stim]    = audioread("SA1.WAV");
fs = Fs_stim;

% TODO: Scale to 65 dB SPL AFTER RIR application (SPL at ear drum not source)
% --> Also take a look at effect of EC

stimdb = 65; % speech level in dB SPL 

HL_freqs = [250 500 1000 2000 4000 6000]; % 6 audiometric frequencies used by HASPI
HL_NH = [0   0   0  0  0  0]; % Normal Hearing
HL_HI = [35  35  40 50 60 65]; % High-freq hearing loss (IEC 60118-15 Moderate HL, Moderately Sloping Group)

%% Design Hearing Aid Gains

fs_h_ha            = 24000;
nfir_24k           = 140;

HL_freqs_design   = [ 0  (HL_freqs ./ (fs_h_ha/2))  1          ];
HL_HI_design      = [ 0  HL_HI                      HL_HI(end) ];
HL_HI_design      = 10 .^ (HL_HI_design ./ 20);
h_ha_mirror_24k   = fir2(nfir_24k, HL_freqs_design, HL_HI_design);

[h_ha_nalr_24k,~] = eb_NALR(HL_HI, nfir_24k, fs_h_ha); %Design the NAL-R filter

[H_ha_nalr_24k, freqs_freqz] = freqz(h_ha_nalr_24k, 1, nfir_24k, fs_h_ha);
[H_ha_mirror_24k,           ~] = freqz(h_ha_mirror_24k, 1, nfir_24k, fs_h_ha);
figure()
plot(HL_freqs, HL_HI, ':k');
hold on;
plot(freqs_freqz, 20*log10(abs(H_ha_nalr_24k)));
plot(freqs_freqz, 20*log10(abs(H_ha_mirror_24k)));
xlabel('Frequency [Hz]')
ylabel('dB')
legend('True Audiogram Mirror', 'NAL-R Filter', 'Mirror Audiogram Design')
title('Hearing Aid Gain Designs')

%% Select Hearing Loss (HL) and Hearing Aid Gain (h_ha)

HL      = HL_NH;
HL_memo = sprintf("HL = [%.0f %.0f %.0f %.0f %.0f %.0f]", HL(1), HL(2), HL(3), HL(4), HL(5), HL(6));

h_ha_list      = {1,                 h_ha_nalr_24k, h_ha_mirror_24k};
h_ha_memo_list = {"No Compensation", "NAL-R",       "Audiogram Mirror"};

h_ha_num  = 1;
h_ha      = h_ha_list{h_ha_num};
h_ha_memo = h_ha_memo_list{h_ha_num};

HA_title_line = sprintf("%s, HA Gain = %s", HL_memo, h_ha_memo);

%% Create results folder

% Results can be saved in the Results Folder
results_dir = sprintf("Results_EC_Test");
if ~exist(results_dir,'Dir')
    mkdir(results_dir)
end

%% Load real RIR data

rir_database = load('../RIR_Databases/rir_databases_4micCombined.mat');
rir_database = rir_database.rir_databases_4micCombined;

if rir_database.fs ~= fs
    error("RIR Database sample rate mismatch");
end

real_RIR_num   = 4;
RIR_list_0deg  = rir_database.rir_data.rir_list_0deg;
RIR_list_90deg = rir_database.rir_data.rir_list_90deg;
G_source_0deg  = RIR_list_0deg{real_RIR_num};
G_source_90deg = RIR_list_90deg{real_RIR_num};
tof_list       = rir_database.rir_data.tof_est{real_RIR_num};
T60_real_RIR   = rir_database.rir_data.T60_msec_list(real_RIR_num) ./ 1000;

% Extract early reflections from real RIR
% er_length      = round((50 / 1000) * Fs_stim); % ~50 msec is "Early Reflections"
% g_rir_left_er  = G_source(:,1);
% g_rir_right_er = G_source(:,2);
% tof_left       = tof_list(1);
% tof_right      = tof_list(2);
% g_rir_left_er  = g_rir_left_er(1:(tof_left   + er_length));
% g_rir_right_er = g_rir_right_er(1:(tof_right + er_length));
% max_er_length  = max([length(g_rir_left_er) length(g_rir_right_er)]);

g_rir_left_0deg  = G_source_0deg(:,1);
g_rir_right_0deg = G_source_0deg(:,2);

g_rir_left_90deg  = G_source_90deg(:,1);
g_rir_right_90deg = G_source_90deg(:,2);


interaural_d  = 15 / 100; % m = cm / 100 

% Spatial noise
S_noise = rir_database.noise_data.noise_list{1}; % Ventillation office room
spatial_room_noise_left  = S_noise(:,1);
spatial_room_noise_right = S_noise(:,2);

%% Generate random RIRs for use as well

% Exponential decay curve
T60_rand_RIR = T60_real_RIR;
N60_rand_RIR = T60_rand_RIR * fs;
L_channel = round(N60_rand_RIR*1.5);
tau = N60_rand_RIR / log(10^3); % -ln(10^(-60/20)) = ln(10^3)
exp_decay = exp(-1 .* (0:(L_channel-1))' ./ tau);

% Option 1: Generate random RIR
g_randn_noise_left   = randn(L_channel,1);
g_randn_noise_right  = randn(L_channel,1);
g_randn_speech_left  = randn(L_channel,1);
g_randn_speech_right = randn(L_channel,1);

% Apply synthetic exponential decay to RIRs
g_randn_noise_left   = g_randn_noise_left(1:L_channel) .* exp_decay;
g_randn_noise_right  = g_randn_noise_right(1:L_channel) .* exp_decay;
g_randn_speech_left  = g_randn_speech_left(1:L_channel) .* exp_decay;
g_randn_speech_right = g_randn_speech_right(1:L_channel) .* exp_decay;

%% NOISE DIRECTION TEST

if 1

SNR_dB     = -12;
DOA_list   = -90:10:90;
HASPI_pre_list  = zeros(length(DOA_list), 1);
HASPI_post_list = zeros(length(DOA_list), 1);

% Generate Uncorrelated noise (Generate outside loop to remove stochasticity
source_noise = randn(length(sref_precalib), 1);

parfor DOA_num = 1:length(DOA_list)

    DOA_noise = DOA_list(DOA_num);
    %fprintf("\n====================================\n")
    fprintf("Noise DOA = %.0f degrees\n", DOA_noise)
    %fprintf("====================================\n")

    % Mic Signals (anechoic)
    DOA_speech = 0; % Front
    [teststim_left_p, teststim_right_p] = HRTF(sref_precalib, DOA_speech, interaural_d, fs);

    % Calibrate dB SPL
    refstim          = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Spatialize Noise
    [noise_left, noise_right] = HRTF(source_noise, DOA_noise, interaural_d, fs);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);
    
    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(DOA_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p, fs, HL, stimdb);
    [HASPI_post_list(DOA_num), ~] = HASPI_v2(refstim, fs, ec_out_p,         fs, HL, stimdb);

end

figure()
plot(DOA_list, HASPI_pre_list, '--b');
hold on
plot(DOA_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Noise azimuth direction [degrees]')
title('Impact of EC on HASPI', {sprintf("Spatial Noise (Synthetic White Noise), T60 = Anechoic at 0 deg, SNR = %.0f dB", SNR_dB), HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_v_DOA.fig', results_dir));

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    SECTION 1: Anechoic Speech                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Anechoic Speech (0 deg) + Anechoic Noise (90 deg)

SNR_dB_list = -12:3:12;
HASPI_pre_list   = zeros(length(SNR_dB_list), 1);
HASPI_post_list  = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (anechoic)
    DOA = 0; % Front
    [teststim_left_p, teststim_right_p] = HRTF(sref_precalib, DOA, interaural_d, fs);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Generate Uncorrelated noise
    source_noise = randn(length(teststim_right_p), 1);

    % Spatialize Noise
     DOA = 90; % Right side
    [noise_left, noise_right] = HRTF(source_noise, DOA, interaural_d, fs);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {"Anechoic Speech (0 deg)", "Anechoic Noise (90 deg white noise)", HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_AnechoicSpeech_AnechoicNoise.fig', results_dir));


%% Anechoic Speech (0 deg) + Reverberant Noise (Real RIR 90 deg)

SNR_dB_list = -12:3:12;
HASPI_pre_list   = zeros(length(SNR_dB_list), 1);
HASPI_post_list  = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (anechoic)
    DOA = 0; % Front
    [teststim_left_p, teststim_right_p] = HRTF(sref_precalib, DOA, interaural_d, fs);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Generate Uncorrelated noise
    source_noise = randn(length(teststim_right_p), 1);

    % Spatialize Noise (Real RIRs)
    noise_left  = filter(g_rir_left_90deg,  1, source_noise);
    noise_right = filter(g_rir_right_90deg, 1, source_noise);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {"Anechoic Speech (0 deg)", sprintf("Reverberant Noise (90 deg white noise, real RIR w/ T60 = %.0f msec)", T60_real_RIR * 1000), HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_AnechoicSpeech_ReverbNoise.fig', results_dir));


%% Anechoic Speech (0 deg) + Spatial Noise Recording


SNR_dB_list   = -12:3:12;
HASPI_pre_list     = zeros(length(SNR_dB_list), 1);
HASPI_post_list     = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (anechoic)
    DOA = 0; % Front
    [teststim_left_p, teststim_right_p] = HRTF(sref_precalib, DOA, interaural_d, fs);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Real Spatial Noise
    noise_right = spatial_room_noise_right(1:length(teststim_right_p));
    noise_left  = spatial_room_noise_left(1:length(teststim_left_p));

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {"Anechoic Speech (0 deg)", "Spatial Noise (Ventilation Noise Recording)", HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_AnechoicSpeech_SpatialNoiseRecording.fig', results_dir));


%% Anechoic Speech (0 deg) + Diffuse Noise


SNR_dB_list   = -12:3:12;
HASPI_pre_list     = zeros(length(SNR_dB_list), 1);
HASPI_post_list     = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (anechoic)
    DOA = 0; % Front
    [teststim_left_p, teststim_right_p] = HRTF(sref_precalib, DOA, interaural_d, fs);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % White noise randomly spatialized to diffuse using random RIRs
    source_noise = randn(length(teststim_right_p), 1);
    noise_left   = filter(g_randn_noise_left,  1, source_noise);
    noise_right  = filter(g_randn_noise_right, 1, source_noise);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);
end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {"Anechoic Speech (0 deg)", sprintf("Diffuse Noise (White noise, Random RIRs, T60 = %.0f msec)", T60_rand_RIR * 1000), HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_AnechoicSpeech_DiffuseNoise.fig', results_dir));


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    SECTION 2: Reverberant Speech                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Reverberant Speech (0 deg) + Anechoic Noise (90 deg)

SNR_dB_list = -12:3:12;
HASPI_pre_list   = zeros(length(SNR_dB_list), 1);
HASPI_post_list  = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (Reverberant, real RIRs)
    teststim_left_p  = filter(g_rir_left_0deg,  1, sref_precalib);
    teststim_right_p = filter(g_rir_right_0deg, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Generate Uncorrelated noise
    source_noise = randn(length(teststim_right_p), 1);

    % Spatialize Noise
     DOA = 90; % Right side
    [noise_left, noise_right] = HRTF(source_noise, DOA, interaural_d, fs);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs); 

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Reverberant Speech (0 deg, Real RIRs T60 = %.0f msec)", T60_real_RIR * 1000), "Anechoic Noise (90 deg white noise)", HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_ReverbSpeech_AnechoicNoise.fig', results_dir));


%% Reverberant Speech (0 deg) + Reverberant Noise (Real RIR 90 deg)

SNR_dB_list = -12:3:12;
HASPI_pre_list   = zeros(length(SNR_dB_list), 1);
HASPI_post_list  = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (Reverberant, real RIRs)
    teststim_left_p  = filter(g_rir_left_0deg,  1, sref_precalib);
    teststim_right_p = filter(g_rir_right_0deg, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Generate Uncorrelated noise
    source_noise = randn(length(teststim_right_p), 1);

    % Spatialize Noise (Real RIRs)
    noise_left  = filter(g_rir_left_90deg,  1, source_noise);
    noise_right = filter(g_rir_right_90deg, 1, source_noise);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);

    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Reverberant Speech (0 deg, Real RIRs T60 = %.0f msec)", T60_real_RIR * 1000), sprintf("Reverberant Noise (90 deg white noise, real RIR w/ T60 = %.0f msec)", T60_real_RIR * 1000), HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_ReverbSpeech_ReverbNoise.fig', results_dir));


%% Reverberant Speech (0 deg) + Spatial Noise Recording


SNR_dB_list   = -12:3:12;
HASPI_pre_list     = zeros(length(SNR_dB_list), 1);
HASPI_post_list     = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (Reverberant, real RIRs)
    teststim_left_p  = filter(g_rir_left_0deg,  1, sref_precalib);
    teststim_right_p = filter(g_rir_right_0deg, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Real Spatial Noise
    noise_right = spatial_room_noise_right(1:length(teststim_right_p));
    noise_left  = spatial_room_noise_left(1:length(teststim_left_p));

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);

    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Reverberant Speech (0 deg, Real RIRs T60 = %.0f msec)", T60_real_RIR * 1000), "Spatial Noise (Ventilation Noise Recording)", HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_ReverbSpeech_SpatialNoiseRecording.fig', results_dir));


%% Reverberant Speech (0 deg) + Diffuse Noise


SNR_dB_list   = -12:3:12;
HASPI_pre_list     = zeros(length(SNR_dB_list), 1);
HASPI_post_list     = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (Reverberant, real RIRs)
    teststim_left_p  = filter(g_rir_left_0deg,  1, sref_precalib);
    teststim_right_p = filter(g_rir_right_0deg, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % White noise randomly spatialized to diffuse using random RIRs
    source_noise = randn(length(teststim_right_p), 1);
    noise_left   = filter(g_randn_noise_left,  1, source_noise);
    noise_right  = filter(g_randn_noise_right, 1, source_noise);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs); 

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);
end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Reverberant Speech (0 deg, Real RIRs T60 = %.0f msec)", T60_real_RIR * 1000), sprintf("Diffuse Noise (White noise, Random RIRs, T60 = %.0f msec)", T60_rand_RIR * 1000), HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_ReverbSpeech_DiffuseNoise.fig', results_dir));



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SECTION 3: Diffuse Speech                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Diffuse Speech (0 deg) + Anechoic Noise (90 deg)

SNR_dB_list = -12:3:12;
HASPI_pre_list   = zeros(length(SNR_dB_list), 1);
HASPI_post_list  = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (randomly spatialized using random RIRs)
    teststim_left_p  = filter(g_randn_speech_left,  1, sref_precalib);
    teststim_right_p = filter(g_randn_speech_right, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Generate Uncorrelated noise
    source_noise = randn(length(teststim_right_p), 1);

    % Spatialize Noise
     DOA = 90; % Right side
    [noise_left, noise_right] = HRTF(source_noise, DOA, interaural_d, fs);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);

    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Diffuse Speech (Random RIRs, T60 =  %.0f msec)", T60_real_RIR * 1000), "Anechoic Noise (90 deg white noise)", HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_DiffuseSpeech_AnechoicNoise.fig', results_dir));



%% Diffuse Speech (0 deg) + Reverberant Noise (Real RIR 90 deg)

SNR_dB_list = -12:3:12;
HASPI_pre_list   = zeros(length(SNR_dB_list), 1);
HASPI_post_list  = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (randomly spatialized using random RIRs)
    teststim_left_p  = filter(g_randn_speech_left,  1, sref_precalib);
    teststim_right_p = filter(g_randn_speech_right, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Generate Uncorrelated noise
    source_noise = randn(length(teststim_right_p), 1);

    % Spatialize Noise (Real RIRs)
    noise_left  = filter(g_rir_left_90deg,  1, source_noise);
    noise_right = filter(g_rir_right_90deg, 1, source_noise);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Diffuse Speech (Random RIRs, T60 = %.0f msec)", T60_real_RIR * 1000), sprintf("Reverberant Noise (90 deg white noise, real RIR w/ T60 = %.0f msec)", T60_real_RIR * 1000), HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_DiffuseSpeech_ReverbNoise.fig', results_dir));


%% Diffuse Speech (0 deg) + Spatial Noise Recording

SNR_dB_list   = -12:3:12;
HASPI_pre_list     = zeros(length(SNR_dB_list), 1);
HASPI_post_list     = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (randomly spatialized using random RIRs)
    teststim_left_p  = filter(g_randn_speech_left,  1, sref_precalib);
    teststim_right_p = filter(g_randn_speech_right, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Real Spatial Noise
    noise_right = spatial_room_noise_right(1:length(teststim_right_p));
    noise_left  = spatial_room_noise_left(1:length(teststim_left_p));

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs); 

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Diffuse Speech (Random RIRs, T60 = %.0f msec)", T60_real_RIR * 1000), "Spatial Noise (Ventilation Noise Recording)", HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_DiffuseSpeech_SpatialNoiseRecording.fig', results_dir));


%% Diffuse Speech (0 deg) + Diffuse Noise

SNR_dB_list   = -12:3:12;
HASPI_pre_list     = zeros(length(SNR_dB_list), 1);
HASPI_post_list     = zeros(length(SNR_dB_list), 1);

parfor SNR_num = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(SNR_num);
    %fprintf("\n====================================\n")
    fprintf("SNR = %.0f dB\n", SNR_dB)
    %fprintf("====================================\n")

    % Mic Signals (randomly spatialized using random RIRs)
    teststim_left_p  = filter(g_randn_speech_left,  1, sref_precalib);
    teststim_right_p = filter(g_randn_speech_right, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % White noise randomly spatialized to diffuse using random RIRs
    source_noise = randn(length(teststim_right_p), 1);
    noise_left   = filter(g_randn_noise_left,  1, source_noise);
    noise_right  = filter(g_randn_noise_right, 1, source_noise);

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs); 

    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(SNR_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(SNR_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);
end

figure()
plot(SNR_dB_list, HASPI_pre_list, '--b');
hold on
plot(SNR_dB_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Test SNR [dB]')
title('Impact of EC on HASPI', {sprintf("Diffuse Speech (Random RIRs, T60 = %.0f msec)", T60_real_RIR * 1000), sprintf("Diffuse Noise (White noise, Random RIRs, T60 = %.0f msec)", T60_rand_RIR * 1000), HA_title_line})
grid on;

saveas(gcf, sprintf('%s/EC_HASPI_DiffuseSpeech_DiffuseNoise.fig', results_dir));


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    SECTION 4: Synthetic T60 TESTS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% T60 Sweep, all synthetic RIRs (random = Diffuse)

SNR_dB_list = [(-12:6:12)' ; 300];% Add 300 for noise free case
T60_list  = [5 10 15 20 25 30 35 40 50 75 100 250 500 750 1000 1250 1500 1750 2000 2250 2500] ./ 1000;

SNR_dB_list = [300]; % temp
T60_list = [5 15 25:25:2500] ./ 1000;%temp
%T60_list = [0.9*1000/fs 5:5:100] ./ 1000;%temp

for SNR_num = 1:length(SNR_dB_list)

    SNR_dB = SNR_dB_list(SNR_num);

    HASPI_pre_list     = zeros(length(T60_list), 1);
    HASPI_post_list     = zeros(length(T60_list), 1);
    
    for T60_num = 1:length(T60_list)
        T60 = T60_list(T60_num);
        %fprintf("\n====================================\n")
        fprintf("T60 = %.0f msec\n", T60*1000)
        %fprintf("====================================\n")
    
        % Exponential decay curve
        N60 = T60 * fs;
        L_channel = round(N60*1.5);
        tau = N60 / log(10^3); % -ln(10^(-30/20)) = ln(10^1.5)
        exp_decay = exp(-1 .* (0:(L_channel-1))' ./ tau);
        
        % Option 1: Generate random RIR
        g_randn_left    = randn(L_channel,1);
        g_randn_right   = randn(L_channel,1);
        g_randn_left_2  = randn(L_channel,1);
        g_randn_right_2 = randn(L_channel,1);
        
        % Apply synthetic exponential decay to RIRs
        g_randn_left    = g_randn_left(1:L_channel) .* exp_decay;
        g_randn_right   = g_randn_right(1:L_channel) .* exp_decay;
        g_randn_left_2  = g_randn_left_2(1:L_channel) .* exp_decay;
        g_randn_right_2 = g_randn_right_2(1:L_channel) .* exp_decay;    
    
        % Mic Signals (randomly spatialized using random RIRs)
        teststim_left_p  = filter(g_randn_left_2,  1, sref_precalib);
        teststim_right_p = filter(g_randn_right_2, 1, sref_precalib);
    
        % Calibrate dB SPL
        refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
        teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
        teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

        % Real Spatial Noise
        noise_right = spatial_room_noise_right(1:length(teststim_right_p));
        noise_left  = spatial_room_noise_left(1:length(teststim_left_p));
    
        % Calibrate SNR
        SNR         = 10 ^ (SNR_dB / 20);
        SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
        noise_right = (noise_right ./ (SNR / SNR_i));
        noise_left  = (noise_left  ./ (SNR / SNR_i));
    
        % Generate signal+noise and signal-noise to separate them later
        teststim_right_p = teststim_right_p + noise_right;
        teststim_left_p  = teststim_left_p  + noise_left;
    
        % Run EC Front-End for signal+noise
        ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  

        % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
        teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
        ec_out_p_ha         = resample(ec_out_p,         24000, fs);
        teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
        ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
        teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
        ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
        
        % Compute HASPI before and after EC
        [HASPI_pre_list(T60_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
        [HASPI_post_list(T60_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

    end
    
    figure()
    plot(T60_list, HASPI_pre_list, '--b');
    hold on
    plot(T60_list, HASPI_post_list, '-b');
    ylim([0 1])
    legend('Before EC', 'After EC')
    ylabel('HASPI')
    xlabel('T60 [msec]')
   % title('Impact of EC on HASPI', {"Diffuse Speech (Random RIRs, Significant 0 deg)", sprintf("Diffuse Noise (White noise, Random RIRs, SNR = %.0f dB)", SNR_dB)}) 
    title('Impact of EC on HASPI', {"Diffuse Speech (Random RIRs)", sprintf("Spatial Noise (Ventilation Noise Recording, SNR = %.0f dB)", SNR_dB)})%, HA_title_line})
    grid on;
    
    saveas(gcf, sprintf('%s/EC_HASPI_v_T60_DiffuseSpeech_DiffuseNoise_SNR_%.0fdB.fig', results_dir, SNR_dB));

end

%% T60 Sweep, all synthetic RIRs (Diffuse but with significant 0 deg energy)

SNR_dB_list = [(-12:6:12)' ; 300];% Add 300 for noise free case
SNR_dB_list = [300]; % temp

T60_list  = [5 10 15 20 25 30 35 40 50 75 100 250 500 750 1000 1250 1500 1750 2000 2250 2500] ./ 1000;

for SNR_num = 1:length(SNR_dB_list)

    SNR_dB = SNR_dB_list(SNR_num);

    HASPI_pre_list     = zeros(length(T60_list), 1);
    HASPI_post_list     = zeros(length(T60_list), 1);
    
    parfor T60_num = 1:length(T60_list)
        T60 = T60_list(T60_num);
        %fprintf("\n====================================\n")
        fprintf("T60 = %.0f msec\n", T60*1000)
        %fprintf("====================================\n")
    
        % Exponential decay curve
        N60 = T60 * fs;
        L_channel = round(N60*1.5);
        tau = N60 / log(10^3); % -ln(10^(-30/20)) = ln(10^1.5)
        exp_decay = exp(-1 .* (0:(L_channel-1))' ./ tau);
        
        % Option 1: Generate random RIR
        g_randn_left    = randn(L_channel,1);
        g_randn_right   = randn(L_channel,1);
        g_randn_left_2  = randn(L_channel,1);
        g_randn_right_2 = randn(L_channel,1);
        
        % Apply synthetic exponential decay to RIRs
        g_randn_left    = g_randn_left(1:L_channel) .* exp_decay;
        g_randn_right   = g_randn_right(1:L_channel) .* exp_decay;
        g_randn_left_2  = g_randn_left_2(1:L_channel) .* exp_decay;
        g_randn_right_2 = g_randn_right_2(1:L_channel) .* exp_decay;  

        % make significant direct sound
        g_randn_right_2(1) = max(abs(g_randn_right_2)) * 4;
        g_randn_left_2(2) = max(abs(g_randn_left_2)) * 4;
    
        % Mic Signals (randomly spatialized using random RIRs)
        teststim_left_p  = filter(g_randn_left_2,  1, sref_precalib);
        teststim_right_p = filter(g_randn_right_2, 1, sref_precalib);
    
        % Calibrate dB SPL
        refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
        teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
        teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

        % Real Spatial Noise
        noise_right = spatial_room_noise_right(1:length(teststim_right_p));
        noise_left  = spatial_room_noise_left(1:length(teststim_left_p));
    
        % Calibrate SNR
        SNR         = 10 ^ (SNR_dB / 20);
        SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
        noise_right = (noise_right ./ (SNR / SNR_i));
        noise_left  = (noise_left  ./ (SNR / SNR_i));
    
        % Generate signal+noise and signal-noise to separate them later
        teststim_right_p = teststim_right_p + noise_right;
        teststim_left_p  = teststim_left_p  + noise_left;
    
        % Run EC Front-End for signal+noise
        ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  
        
        % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
        teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
        ec_out_p_ha         = resample(ec_out_p,         24000, fs);
        teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
        ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
        teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
        ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
        
        % Compute HASPI before and after EC
        [HASPI_pre_list(T60_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
        [HASPI_post_list(T60_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

    end
    
    figure()
    plot(T60_list, HASPI_pre_list, '--b');
    hold on
    plot(T60_list, HASPI_post_list, '-b');
    ylim([0 1])
    legend('Before EC', 'After EC')
    ylabel('HASPI')
    xlabel('T60 [msec]')
   % title('Impact of EC on HASPI', {"Diffuse Speech (Random RIRs, Significant 0 deg)", sprintf("Diffuse Noise (White noise, Random RIRs, SNR = %.0f dB)", SNR_dB)}) 
    title('Impact of EC on HASPI', {"Diffuse Speech (Random RIRs, Significant 0 deg)", sprintf("Spatial Noise (Ventilation Noise Recording, SNR = %.0f dB)", SNR_dB), HA_title_line})
    grid on;
    
    saveas(gcf, sprintf('%s/EC_HASPI_v_T60_DiffuseSpeech_DiffuseNoise_SNR_%.0fdB.fig', results_dir, SNR_dB));

end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    SECTION 5: T60 TESTS  w/ Real RIRs                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Load Real RIR Data

rir_database = load('../RIR_Databases/rir_databases_4micCombined.mat');
rir_database = rir_database.rir_databases_4micCombined;

if rir_database.fs ~= Fs_stim
    error("RIR Database sample rate mismatch");
end

RIR_list_0deg  = rir_database.rir_data.rir_list_0deg;
RIR_list_90deg = rir_database.rir_data.rir_list_90deg;
room_memos     = rir_database.rir_data.rir_memos;
T60_msec_list  = rir_database.rir_data.T60_msec_list;
T60_list       = T60_msec_list ./ 1000;
tof_list       = rir_database.rir_data.tof_est;

% Reduce to subset
%T60_list = T60_list(1);



%% T60 Sweep, all Real RIRs


SNR_dB_list = [(-12:6:12)' ; 300];% Add 300 for noise free case
SNR_dB_list = [300]; % temp

for SNR_num = 1:length(SNR_dB_list)

    SNR_dB = SNR_dB_list(SNR_num);

    HASPI_pre_list     = zeros(length(T60_list), 1);
    HASPI_post_list     = zeros(length(T60_list), 1);
    
    parfor T60_num = 1:length(T60_list)
        T60 = T60_list(T60_num);
        %fprintf("\n====================================\n")
        fprintf("T60 = %.0f msec\n", T60*1000)
        %fprintf("====================================\n")

        % Real RIRs
        G_source        = RIR_list_0deg{T60_num};
        g_left_2  = G_source(:,1);
        g_right_2 = G_source(:,2);
    
        % Mic Signals (randomly spatialized using random RIRs)
        teststim_left_p  = filter(g_left_2,  1, sref_precalib);
        teststim_right_p = filter(g_right_2, 1, sref_precalib);
    
        % Calibrate dB SPL
        refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
        teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
        teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);
    
        % Real Spatial Noise
        noise_right = spatial_room_noise_right(1:length(teststim_right_p));
        noise_left  = spatial_room_noise_left(1:length(teststim_left_p));
    
        % Calibrate SNR
        SNR         = 10 ^ (SNR_dB / 20);
        SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
        noise_right = (noise_right ./ (SNR / SNR_i));
        noise_left  = (noise_left  ./ (SNR / SNR_i));
    
        % Generate signal+noise and signal-noise to separate them later
        teststim_right_p = teststim_right_p + noise_right;
        teststim_left_p  = teststim_left_p  + noise_left;
    
        % Run EC Front-End for signal+noise
        ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  
        
        % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
        teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
        ec_out_p_ha         = resample(ec_out_p,         24000, fs);
        teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
        ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
        teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
        ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
        
        % Compute HASPI before and after EC
        [HASPI_pre_list(T60_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
        [HASPI_post_list(T60_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

    end
    
    figure()
    plot(T60_list, HASPI_pre_list, '--ob');
    hold on
    plot(T60_list, HASPI_post_list, '-ob');
    ylim([0 1])
    legend('Before EC', 'After EC')
    ylabel('HASPI')
    xlabel('T60 [msec]')
    title('Impact of EC on HASPI', {"Reverberant Speech (Real RIRs)", sprintf("Spatial Noise (Ventilation Noise Recording, SNR = %.0f dB)", SNR_dB), HA_title_line})
    grid on;
    
    saveas(gcf, sprintf('%s/EC_HASPI_v_T60_ReverberantSpeech_DiffuseNoise_SNR_%.0fdB.fig', results_dir, SNR_dB));

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         SECTION 6: T60 TESTS  w/ truncated MYRiAD 2.1 sec RIR           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% Real RIRs
room_num     = 4;
G_source     = RIR_list_0deg{room_num};
g_left_SAL  = G_source(:,1);
g_right_SAL = G_source(:,2);

SNR_dB_list = [(-12:6:12)' ; 300];% Add 300 for noise free case
T60_list = [5 15 0:100:2100] ./ 1000;

SNR_dB_list = [300]; % Temp

for SNR_num = 1:length(SNR_dB_list)

    SNR_dB = SNR_dB_list(SNR_num);

    HASPI_pre_list     = zeros(length(T60_list), 1);
    HASPI_post_list     = zeros(length(T60_list), 1);
    
    parfor T60_num = 1:length(T60_list)
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
        teststim_left_p  = filter(g_left,  1, sref_precalib);
        teststim_right_p = filter(g_right, 1, sref_precalib);
    
        % Calibrate dB SPL
        teststim_calib_gain = 1 / max([rms(teststim_left_p) rms(teststim_right_p)]) * 20e-6*10^(stimdb/20);
        refstim        = sref_precalib  / rms(sref_precalib)  * 20e-6*10^(stimdb/20);
        teststim_left_p  = teststim_left_p  * teststim_calib_gain;
        teststim_right_p = teststim_right_p * teststim_calib_gain;
    
        % Add noise
        SNR = 10 ^ (SNR_dB / 20);
        s_noise_left  = randn(length(teststim_left_p),   1);
        s_noise_right = randn(length(teststim_right_p),   1);
        SNR_i       = max([rms(teststim_left_p) rms(teststim_right_p)]) / rms(s_noise_left);
        teststim_left_p   = teststim_left_p  + (s_noise_left ./ (SNR / SNR_i));
        teststim_right_p  = teststim_right_p + (s_noise_right ./ (SNR / SNR_i));
    
        % Run EC Front-End for signal+noise
        ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  
        
        % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
        % teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
        % ec_out_p_ha         = resample(ec_out_p,         24000, fs);
        % teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
        % ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
        % teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
        % ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);

        % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
        [ec_out_p_ha, fsamp] = eb_Resamp24kHz(ec_out_p,  Fs_stim);
        nsamp       = length(ec_out_p_ha);
        nfir        = length(h_ha) - 1;
        ec_out_p_ha = conv(ec_out_p_ha,  h_ha);
        ec_out_p_ha = ec_out_p_ha(nfir+1:nfir+nsamp);
        ec_out_p_ha = resample(ec_out_p_ha,  Fs_stim, fsamp);


        [teststim_right_p_ha, fsamp] = eb_Resamp24kHz(teststim_left_p,  Fs_stim);
        nsamp               = length(teststim_right_p_ha);
        nfir                = length(h_ha) - 1;
        teststim_right_p_ha = conv(teststim_right_p_ha,  h_ha);
        teststim_right_p_ha = teststim_right_p_ha(nfir+1:nfir+nsamp);
        teststim_right_p_ha = resample(teststim_right_p_ha,  Fs_stim, fsamp);
        
        % Compute HASPI before and after EC
        [HASPI_pre_list(T60_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
        [HASPI_post_list(T60_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

    end
    
    figure()
    plot(T60_list, HASPI_pre_list, '--b');
    hold on
    plot(T60_list, HASPI_post_list, '-b');
    ylim([0 1])
    legend('Before EC', 'After EC')
    ylabel('HASPI')
    xlabel('T60 [msec]')
   % title('Impact of EC on HASPI', {"Diffuse Speech (Random RIRs, Significant 0 deg)", sprintf("Diffuse Noise (White noise, Random RIRs, SNR = %.0f dB)", SNR_dB)}) 
    title('Impact of EC on HASPI', {"Reverberant Speech (MYRiAD 2.1 sec RIR truncated exponentially)", sprintf("Spatial Noise (Ventilation Noise Recording, SNR = %.0f dB)", SNR_dB), HA_title_line})
    grid on;
    
    saveas(gcf, sprintf('%s/EC_HASPI_v_T60_DiffuseSpeech_DiffuseNoise_SNR_%.0fdB.fig', results_dir, SNR_dB));

end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         SECTION 8: T30 = 2sec with varying levels of direct sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
direct_gain_db_list = 0:1:42; % dB


SNR_dB = 300; % noise free
T60 = 2;


HASPI_pre_list  = zeros(length(direct_gain_db_list), 1);
HASPI_post_list = zeros(length(direct_gain_db_list), 1);

parfor direct_gain_num = 1:length(direct_gain_db_list)
    %fprintf("\n====================================\n")
    fprintf("T60 = %.0f msec\n", T60*1000)
    %fprintf("====================================\n")

    direct_gain = 10^(direct_gain_db_list(direct_gain_num)/20);

    % Exponential decay curve
    N60 = T60 * fs;
    L_channel = round(N60*1.5);
    tau = N60 / log(10^3); % -ln(10^(-30/20)) = ln(10^1.5)
    exp_decay = exp(-1 .* (0:(L_channel-1))' ./ tau);
    
    % Option 1: Generate random RIR
    g_randn_left    = randn(L_channel,1);
    g_randn_right   = randn(L_channel,1);
    g_randn_left_2  = randn(L_channel,1);
    g_randn_right_2 = randn(L_channel,1);
    
    % Apply synthetic exponential decay to RIRs
    g_randn_left    = g_randn_left(1:L_channel) .* exp_decay;
    g_randn_right   = g_randn_right(1:L_channel) .* exp_decay;
    g_randn_left_2  = g_randn_left_2(1:L_channel) .* exp_decay;
    g_randn_right_2 = g_randn_right_2(1:L_channel) .* exp_decay;  

    % make significant direct sound
    g_randn_right_2(1) = max(abs(g_randn_right_2)) * direct_gain;
    g_randn_left_2(2) = max(abs(g_randn_left_2)) * direct_gain;

    addpath ../dereverb/utilities/matlab/
    figure()
    edc = EDC(g_randn_left_2);
    subplot(2,1,1)
    plot((0:(length(g_randn_left_2)-1)) .* (1/fs), g_randn_left_2);
    grid on;
    title('IR')
    xlabel('sec')
    xlim([-0.05 2])
    subplot(2,1,2)
    plot((0:(length(edc)-1)) .* (1/fs), 10*log10(edc));
    grid on;
    title('EDC')
    ylabel('dB')
    xlabel('sec')
    xlim([-0.05 2])
    ylim([-80 6])
    sgtitle(sprintf("Synthetic RIR (T30 = 2 sec, Direct gain = %.0f dB)", 20*log10(direct_gain)))

    % Mic Signals (randomly spatialized using random RIRs)
    teststim_left_p  = filter(g_randn_left_2,  1, sref_precalib);
    teststim_right_p = filter(g_randn_right_2, 1, sref_precalib);

    % Calibrate dB SPL
    refstim        = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
    teststim_left_p  = teststim_left_p/rms(teststim_left_p)*20e-6*10^(stimdb/20);
    teststim_right_p = teststim_right_p/rms(teststim_right_p)*20e-6*10^(stimdb/20);

    % Real Spatial Noise
    noise_right = spatial_room_noise_right(1:length(teststim_right_p));
    noise_left  = spatial_room_noise_left(1:length(teststim_left_p));

    % Calibrate SNR
    SNR         = 10 ^ (SNR_dB / 20);
    SNR_i       = rms(teststim_right_p) / max([rms(noise_left) rms(noise_right)]);
    noise_right = (noise_right ./ (SNR / SNR_i));
    noise_left  = (noise_left  ./ (SNR / SNR_i));

    % Generate signal+noise and signal-noise to separate them later
    teststim_right_p = teststim_right_p + noise_right;
    teststim_left_p  = teststim_left_p  + noise_left;

    % Run EC Front-End for signal+noise
    ec_out_p  = EC_FrontEnd(teststim_left_p, teststim_right_p, fs);  
    
    % Apply hearing aid gain Filter to compensate hearing loss (if HL = zeros, filter will be bypass)
    teststim_right_p_ha = resample(teststim_right_p, 24000, fs);
    ec_out_p_ha         = resample(ec_out_p,         24000, fs);
    teststim_right_p_ha = filter(h_ha, 1, teststim_right_p_ha);
    ec_out_p_ha         = filter(h_ha, 1, ec_out_p_ha);
    teststim_right_p_ha = resample(teststim_right_p_ha, fs, 24000);
    ec_out_p_ha         = resample(ec_out_p_ha,         fs, 24000);
    
    % Compute HASPI before and after EC
    [HASPI_pre_list(direct_gain_num), ~]  = HASPI_v2(refstim, fs, teststim_right_p_ha, fs, HL, stimdb);
    [HASPI_post_list(direct_gain_num), ~] = HASPI_v2(refstim, fs, ec_out_p_ha,         fs, HL, stimdb);

end

figure()
plot(direct_gain_db_list, HASPI_pre_list, '--b');
hold on
plot(direct_gain_db_list, HASPI_post_list, '-b');
ylim([0 1])
legend('Before EC', 'After EC')
ylabel('HASPI')
xlabel('Direct Sound Gain [dB]')
% title('Impact of EC on HASPI', {"Diffuse Speech (Random RIRs, Significant 0 deg)", sprintf("Diffuse Noise (White noise, Random RIRs, SNR = %.0f dB)", SNR_dB)}) 
title('Impact of EC on HASPI', {"Diffuse Speech (Random RIRs, Significant 0 deg)", sprintf("Spatial Noise (Ventilation Noise Recording, SNR = %.0f dB)", SNR_dB), HA_title_line})
grid on;
