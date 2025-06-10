restoredefaultpath

close all
clear
clc

myCluster = parcluster('Processes');
delete(myCluster.Jobs)
delete(gcp('nocreate'))
pause(5)
parpool


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
rng(724); % Reproduce specific seed


%% Init

SNR_dB = 300;

HL_freqs = [250 500 1000 2000 4000 6000]; % 6 audiometric frequencies used by HASPI
HL_NH    = [0   0   0  0  0  0]; % Normal Hearing
HL_HI    = [35  35  40 50 60 65]; % High-freq hearing loss (IEC 60118-15 Moderate HL, Moderately Sloping Group)

[sref_precalib, Fs_stim] = audioread("SA1.WAV");

stimdb = 65; % speech level in dB SPL 

my_colours = get(gca,'colororder');



%% Create results folder

% Results can be saved in the Results Folder
results_dir_parent = sprintf("Results_HA_Gain_Comparison_T60Sweep_RealER_%.0fdB_SNR", SNR_dB);
if ~exist(results_dir_parent,'Dir')
    mkdir(results_dir_parent)
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
plot(freqs_freqz, 20*log10(abs(H_ha_nocomp_24k)), 'Color', my_colours(1,:));
plot(freqs_freqz, 20*log10(abs(H_ha_mirror_24k)), 'Color', my_colours(2,:));
plot(freqs_freqz, 20*log10(abs(H_ha_hybrid_24k)), 'Color', my_colours(3,:));
plot(freqs_freqz, 20*log10(abs(H_ha_nalr_24k)),   'Color', my_colours(4,:));
xlabel('Frequency [Hz]')
ylabel('dB')
legend('True Audiogram Mirror', 'No Compensation', 'Audiogram Mirror Design', 'Hybrid', 'NAL-R Filter')
title('Hearing Aid Gain Designs')

saveas(gcf, sprintf('%s/HearingAidGains.fig', results_dir_parent));

% Create Test vectors
h_ha_list           = { h_ha_nocomp_24k, h_ha_mirror_24k,         h_ha_hybrid_24k, h_ha_nalr_24k   };
h_ha_memos          = { "No Gain",       "Mirror Audiogram Gain", "Hybrid Gain",   "NAL-R Gain"    }; 
h_ha_memos_nospaces = { "No_Gain",       "Mirror_Audiogram_Gain", "Hybrid_Gain",   "NALR_Gain"     }; 
line_colour_list    = { my_colours(1,:), my_colours(2,:),         my_colours(3,:), my_colours(4,:) };
line_marker_list    = { "o",             "^",                     "square",        "diamond"       };

%% Load real RIR data

rir_database = load('../RIR_Databases/rir_databases_4micCombined_rirOnly.mat');
rir_database = rir_database.rir_databases_4micCombined_rir_only;

if rir_database.fs ~= Fs_stim
    error("RIR Database sample rate mismatch");
end

room_num       = 2;
%RIR_list_0deg  = rir_database.rir_data.rir_list_0deg;
G_source       = rir_database.rir_data.rir_list_0deg{room_num};
tof_list       = rir_database.rir_data.tof_est{room_num};

% Extract early reflections from real RIR
er_length      = round((50 / 1000) * Fs_stim); % ~50 msec is "Early Reflections"
g_rir_er       = G_source(:,1);
tof            = tof_list(1);
g_rir_er       = g_rir_er(1:(tof + er_length));
er_full_length = length(g_rir_er);

T60_list = (0:250:2500) ./ 1000;

% Pre-generate the longest T60 and truncate later to reduce variance in results
L_channel_max = round(max(T60_list) * Fs_stim * 2);
g_tail_max = randn(L_channel_max,1);


%% Hearing Aid Gain Loop

for ha_gain_num = 1:length(h_ha_list)

    results_dir = sprintf("%s/%s", results_dir_parent, h_ha_memos_nospaces{ha_gain_num});
    % Results can be saved in the Results Folder
    if ~exist(results_dir,'Dir')
        mkdir(results_dir)
    end

    h_ha = h_ha_list{ha_gain_num};
    
    % T60 Loop
     
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
        
        fs = Fs_stim;
        
        if T60 == 0
            g_rir  = 1;
        else
            
            % Exponential decay curve
            N60               = T60 * fs;
            tau               = N60 / log(10^3); % -ln(10^(-60/20)) = ln(10^3)
            lr_start          = tof + round((11 / 1000) * fs); % Give 15 msec before starting exponential decay
            L_channel         = max(lr_start + round(N60*1.5));
            tail_length_list  = L_channel - lr_start + 1;
            tail_length       = tail_length_list(1);
            exp_decay         = exp(-1 .* (1:tail_length)'  ./ tau);
            
            % Apply synthetic exponential decay to tail
            g_rir_tail  = g_tail_max(1:tail_length)   .* exp_decay;
            g_rir_tail  = [zeros((L_channel - length(g_rir_tail)),  1) ; g_rir_tail];
            g_rir_tail  = g_rir_tail  .* 0.2 * (max(abs(g_rir_er))  / max(abs(g_rir_tail)));
            
            % Zero pad ER to same length 
            if length(g_rir_er) < L_channel
                g_rir_er_zp  = [g_rir_er  ; zeros((L_channel - length(g_rir_er)),  1)];
            else
                g_rir_er_zp  = g_rir_er(1:L_channel);
            end
        
            % Construct full RIR: Add early reflections to tail
            g_rir  = g_rir_er_zp  + g_rir_tail;
    
        end

        % Plot
        % addpath ../dereverb/utilities/matlab/
        % edc = EDC(g_rir);
        % figure()
        % subplot(3,1,1)
        % plot((0:(length(G_source(:,1))-1)) .* (1/fs), G_source(:,1));
        % grid on;
        % subplot(3,1,2)
        % plot((0:(length(g_rir)-1)) .* (1/fs), g_rir);
        % grid on;
        % subplot(3,1,3)
        % plot((0:(length(edc)-1)) .* (1/fs), 10*log10(edc));
        % grid on;
        % sgtitle(sprintf("T60 + %.0f msec", T60 * 1000))
    
        % Mic Signals
        teststim  = filter(g_rir,  1, sref_precalib);
    
        % Calibrate dB SPL
        refstim   = sref_precalib/rms(sref_precalib)*20e-6*10^(stimdb/20);
        teststim  = teststim/rms(teststim)*20e-6*10^(stimdb/20);
    
        % Add noise
        SNR      = 10 ^ (SNR_dB / 20);
        s_noise  = randn(length(teststim),   1);
        SNR_i    = rms(teststim) / rms(s_noise);
        teststim = teststim  + (s_noise ./ (SNR / SNR_i));
    
        % Compute metrics
        [metrics_HI, nalr_struct_HI] = run_all_metrics(refstim, teststim, Fs_stim, HL_HI, stimdb, h_ha);
        
        % Save NAL-R Filter info
        audiowrite(sprintf("%s/test_stim_noNALR_T60_%.0fmsec.wav",   results_dir, T60*1000), teststim .* (0.5 / max(abs(teststim))), Fs_stim);
        audiowrite(sprintf("%s/test_stim_withNALR_T60_%.0fmsec.wav", results_dir, T60*1000), nalr_struct_HI.teststim_post_nalr .* (0.5 / max(abs(nalr_struct_HI.teststim_post_nalr))), Fs_stim);
        save(sprintf("%s/nalr_struct_HI_T60_%.0fmsec.mat", results_dir, T60*1000),"-fromstruct", nalr_struct_HI);

        STMI_list_HI(T60_num)           = metrics_HI.STMI;
        NSIM_FT_list_HI(T60_num)        = metrics_HI.NSIM_FT;
        NSIM_MR_list_HI(T60_num)        = metrics_HI.NSIM_MR;
        HASPI_list_HI(T60_num)          = metrics_HI.HASPI;
        STOI_list_HI(T60_num)           = metrics_HI.STOI;
        HASQI_list_HI(T60_num)          = metrics_HI.HASQI;
        VISQOL_list_HI(T60_num)         = metrics_HI.VISQOL;
        
    end
    
    test_memo_NH = sprintf("(%s, Synthetic RIR w/ Real ERs, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f])", h_ha_memos{ha_gain_num}, SNR_dB, stimdb, HL_NH(1), HL_NH(2), HL_NH(3), HL_NH(4), HL_NH(5), HL_NH(6));
    test_memo_HI = sprintf("(%s, Synthetic RIR w/ Real ERs, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f])", h_ha_memos{ha_gain_num}, SNR_dB, stimdb, HL_HI(1), HL_HI(2), HL_HI(3), HL_HI(4), HL_HI(5), HL_HI(6));
    
    % figure()
    % subplot(1,2,1)
    % plot(T60_list, STMI_list_NH, "-or")
    % hold on;
    % plot(T60_list, NSIM_FT_list_NH, "-^g")
    % plot(T60_list, NSIM_MR_list_NH, "-^c")
    % plot(T60_list, HASPI_list_NH, "-squareb")
    % plot(T60_list, STOI_list_NH, "-diamondm")
    % xlabel('T60 (sec)')
    % legend("STMI", "NSIM FT", "NSIM MR", "HASPI", "STOI")
    % title("Speech Intelligibility Predictors", test_memo_NH)
    % 
    % subplot(1,2,2)
    % plot(T60_list, STMI_list_HI, "-or")
    % hold on;
    % plot(T60_list, NSIM_FT_list_HI, "-^g")
    % plot(T60_list, NSIM_MR_list_HI, "-^c")
    % plot(T60_list, HASPI_list_HI, "-squareb")
    % xlabel('T60 (sec)')
    % legend("STMI", "NSIM FT", "NSIM MR", "HASPI")
    % title("Speech Intelligibility Predictors", test_memo_HI)
    % 
    % saveas(gcf, sprintf('%s/SI_v_T60.fig', results_dir));
    % 
    % figure()
    % subplot(1,2,1)
    % plot(T60_list, HASQI_list_NH, "-squareb")
    % hold on;
    % plot(T60_list, VISQOL_list_NH, "-diamondm")
    % xlabel('T60 (sec)')
    % legend("HASQI", "VISQOL")
    % title("Speech Quality Predictors", test_memo_NH)
    % 
    % subplot(1,2,2)
    % plot(T60_list, HASQI_list_HI, "-squareb")
    % xlabel('T60 (sec)')
    % legend("HASQI")
    % title("Speech Quality Predictors", test_memo_HI)
    % 
    % saveas(gcf, sprintf('%s/SQ_v_T60.fig', results_dir));
    
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

    % Plot MR/FT NSIM, STMI and HASPI for all hearing aid gains on the same plot
    figure(10)
    plot(T60_list, STMI_list_HI, 'Color', line_colour_list{ha_gain_num}, 'Marker', line_marker_list{ha_gain_num});
    hold on;

    figure(11)
    plot(T60_list, NSIM_FT_list_HI, 'Color', line_colour_list{ha_gain_num}, 'Marker', line_marker_list{ha_gain_num});
    hold on;

    figure(12)
    plot(T60_list, NSIM_MR_list_HI, 'Color', line_colour_list{ha_gain_num}, 'Marker', line_marker_list{ha_gain_num});
    hold on;

    figure(13)
    plot(T60_list, HASPI_list_HI, 'Color', line_colour_list{ha_gain_num}, 'Marker', line_marker_list{ha_gain_num});
    hold on;

    figure(14)
    plot(T60_list, HASQI_list_HI, 'Color', line_colour_list{ha_gain_num}, 'Marker', line_marker_list{ha_gain_num});
    hold on;
end

test_memo_NH = sprintf("(Synthetic RIR w/ Real ERs, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f])", SNR_dB, stimdb, HL_NH(1), HL_NH(2), HL_NH(3), HL_NH(4), HL_NH(5), HL_NH(6));
test_memo_HI = sprintf("(Synthetic RIR w/ Real ERs, SNR = %.0f dB, Stimulus = %.0f dBSPL, HL = [%.0f %.0f %.0f %.0f %.0f %.0f])", SNR_dB, stimdb, HL_HI(1), HL_HI(2), HL_HI(3), HL_HI(4), HL_HI(5), HL_HI(6));


figure(10)
xlabel('T60 (sec)')
ylabel('STMI')
legend("No Gain", "Mirror Audiogram Gain", "Hybrid Gain", "NAL-R Gain");
title("STMI for various hearing aid gains", test_memo_HI)

saveas(gcf, sprintf('%s/STMI_v_T60_variousHAGains.fig', results_dir_parent));

figure(11)
xlabel('T60 (sec)')
ylabel('FT NSIM')
legend("No Gain", "Mirror Audiogram Gain", "Hybrid Gain", "NAL-R Gain");
title("FT NSIM for various hearing aid gains", test_memo_HI)

saveas(gcf, sprintf('%s/NSIM_FT_v_T60_variousHAGains.fig', results_dir_parent));

figure(12)
xlabel('T60 (sec)')
ylabel('MR NSIM')
legend("No Gain", "Mirror Audiogram Gain", "Hybrid Gain", "NAL-R Gain");
title("MR NSIM for various hearing aid gains", test_memo_HI)

saveas(gcf, sprintf('%s/NSIM_MR_v_T60_variousHAGains.fig', results_dir_parent));

figure(13)
xlabel('T60 (sec)')
ylabel('HASPI')
legend("No Gain", "Mirror Audiogram Gain", "Hybrid Gain", "NAL-R Gain");
title("HASPI for various hearing aid gains", test_memo_HI)

saveas(gcf, sprintf('%s/HASPI_v_T60_variousHAGains.fig', results_dir_parent));

figure(14)
xlabel('T60 (sec)')
ylabel('HASQI')
legend("No Gain", "Mirror Audiogram Gain", "Hybrid Gain", "NAL-R Gain");
title("HASQI for various hearing aid gains", test_memo_HI)

saveas(gcf, sprintf('%s/HASQI_v_T60_variousHAGains.fig', results_dir_parent));