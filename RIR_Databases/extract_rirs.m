restoredefaultpath

close all
clear
clc


addpath('HRIR_Universitat_Oldenburg/HRIR_database_mat/')
addpath('HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/')
addpath('HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/courtyard/')
addpath('HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/cafeteria/')
addpath('MYRiAD_V2/tools/MATLAB')

addpath('../dereverb/utilities/matlab/')

enable_plots = true;

fig_size = [10 10 1300 600];

%% Params

fs_output = 16000;
rir_databases.fs = fs_output;
rir_databases_4micCombined.fs = fs_output;

noise_num_4mic = 1;
rir_num_4mic   = 1;

%% HRIR Database

rir_databases.hrir.channel_memos = ["BTE Front Left" "BTE Front Right" "BTE Middle Left" "BTE Middle Right" "BTE Rear Left" "BTE Rear Right"];
rir_databases_4micCombined.channel_memos = ["BTE Front Left" "BTE Front Right" "BTE Rear Left" "BTE Rear Right"];

% rooms               = {'Office_II', 'Courtyard',   'Cafeteria' };
% speaker_locs_0deg   = {'A',         'B',           'A'         };
% speaker_locs_90deg  = {'C',         'D',           'C'         };
% hear_orientations   = {2,            1,             1          };
% rir_memos           = {"office",    "courtyard",   "cafeteria"   };
% cited_T60_msec_list = {300,         900,           1250        }; % Paper
% T60_msec_list       = {400,         390,           542         }; % My measurement
% data_set            = 'bte';

rooms               = { 'Courtyard', 'Office_II', 'Cafeteria'  };
speaker_locs_0deg   = { 'B',         'A',          'A'         };
speaker_locs_90deg  = { 'D',         'C',          'C'         };
hear_orientations   = { 1,           2,            1           };
rir_memos           = { "courtyard", "office",     "cafeteria" };
cited_T60_msec_list = { 900,         300,          1250        }; % Paper
T60_msec_list       = { 390,         400,          542         }; % My measurement
data_set            = 'bte';

for rir_num = 1:length(rir_memos)

    room              = rooms{rir_num};
    head_orientation  = hear_orientations{rir_num};
    speaker_loc_0deg  = speaker_locs_0deg{rir_num};
    speaker_loc_90deg = speaker_locs_90deg{rir_num};
    
    HRIR_data_0deg  = loadHRIR(room, head_orientation, speaker_loc_0deg,  data_set);
    HRIR_data_90deg = loadHRIR(room, head_orientation, speaker_loc_90deg, data_set);
    
    if HRIR_data_0deg.fs ~= fs_output
        [resample_p, resample_q] = rat(fs_output / HRIR_data_0deg.fs);
        HRIR_data_0deg.data  = resample(HRIR_data_0deg.data,  resample_p, resample_q);
        HRIR_data_90deg.data = resample(HRIR_data_90deg.data, resample_p, resample_q);
    end

    % Time-align RIRs and
    % Remove leading measurement noise manually (dont want to convolve this, it isnt really part of the RIR)
    [~, M] = size(HRIR_data_0deg.data);
    tofs_est_0deg = zeros(M,1);
    tofs_est_90deg = zeros(M,1);
    for ch = 1:M

        [HRIR_data_0deg.data(:,ch),  tofs_est_0deg(ch)]  = remove_air_meas_noise(HRIR_data_0deg.data(:,ch));
        [HRIR_data_90deg.data(:,ch), tofs_est_90deg(ch)] = remove_air_meas_noise(HRIR_data_90deg.data(:,ch));

        % Compensate delay --> TODO Should only happen in the DAP algo i think (dont wanna get rid of ITDs)
        %tdf = [zeros(tofs_est(ch), 1) ; 1];
        %HRIR_data.data(:,ch) = filter(tdf, 1, HRIR_data.data(:,ch));
    end
    
    % Channels:
    % - 1: Left Front
    % - 2: Right Front
    % - 3: Left Middle
    % - 4: Right Middle
    % - 5: Left Rear
    % - 6: Front Rear
    rir_left_frontmic_0deg   = HRIR_data_0deg.data(:,1);
    rir_right_frontmic_0deg  = HRIR_data_0deg.data(:,2);
    rir_left_middlemic_0deg  = HRIR_data_0deg.data(:,3);
    rir_right_middlemic_0deg = HRIR_data_0deg.data(:,4);
    rir_left_rearmic_0deg    = HRIR_data_0deg.data(:,5);
    rir_right_rearmic_0deg   = HRIR_data_0deg.data(:,6);

    rir_left_frontmic_90deg   = HRIR_data_90deg.data(:,1);
    rir_right_frontmic_90deg  = HRIR_data_90deg.data(:,2);
    rir_left_middlemic_90deg  = HRIR_data_90deg.data(:,3);
    rir_right_middlemic_90deg = HRIR_data_90deg.data(:,4);
    rir_left_rearmic_90deg    = HRIR_data_90deg.data(:,5);
    rir_right_rearmic_90deg   = HRIR_data_90deg.data(:,6);

    rir_databases.hrir.rir_data.rir_list_0deg{rir_num} = [rir_left_frontmic_0deg    ...
                                                          rir_right_frontmic_0deg  ...
                                                          rir_left_middlemic_0deg  ...
                                                          rir_right_middlemic_0deg ...
                                                          rir_left_rearmic_0deg    ...
                                                          rir_right_rearmic_0deg ];

    rir_databases.hrir.rir_data.rir_list_90deg{rir_num} = [rir_left_frontmic_90deg    ...
                                                           rir_right_frontmic_90deg  ...
                                                           rir_left_middlemic_90deg  ...
                                                           rir_right_middlemic_90deg ...
                                                           rir_left_rearmic_90deg    ...
                                                           rir_right_rearmic_90deg ];

    rir_databases.hrir.rir_data.rir_memos(rir_num)     = rir_memos{rir_num};
    rir_databases.hrir.rir_data.T60_msec_list(rir_num) = T60_msec_list{rir_num};
    rir_databases.hrir.rir_data.tof_est_0deg{rir_num}  = tofs_est_0deg;
    rir_databases.hrir.rir_data.tof_est_90deg{rir_num} = tofs_est_90deg;

    % Save 4-mic subset to a combined structures
    rir_databases_4micCombined.rir_data.rir_list_0deg{rir_num} = [rir_left_frontmic_0deg    ...
                                                                       rir_right_frontmic_0deg  ...
                                                                       rir_left_rearmic_0deg    ...
                                                                       rir_right_rearmic_0deg ];

    rir_databases_4micCombined.rir_data.rir_list_90deg{rir_num} = [rir_left_frontmic_90deg    ...
                                                                        rir_right_frontmic_90deg  ...
                                                                        rir_left_rearmic_90deg    ...
                                                                        rir_right_rearmic_90deg ];

    rir_databases_4micCombined.rir_data.rir_memos(rir_num_4mic)     = rir_memos{rir_num};
    rir_databases_4micCombined.rir_data.T60_msec_list(rir_num_4mic) = T60_msec_list{rir_num};
    rir_databases_4micCombined.rir_data.tof_est_0deg{rir_num_4mic}  = tofs_est_0deg([1 2 5 6]);
    rir_databases_4micCombined.rir_data.tof_est_90deg{rir_num_4mic} = tofs_est_90deg([1 2 5 6]);

    if enable_plots
        edc_1 = EDC(rir_left_frontmic_0deg);
        edc_2 = EDC(rir_right_frontmic_0deg);
        edc_3 = EDC(rir_left_middlemic_0deg);
        edc_4 = EDC(rir_right_middlemic_0deg);
        edc_5 = EDC(rir_left_rearmic_0deg);
        edc_6 = EDC(rir_right_rearmic_0deg);
    
        figure()
        set(gcf,'Position',fig_size)
        subplot(1,2,1)
        plot((0:(length(rir_left_frontmic_0deg)-1)) .* (1/fs_output), rir_left_frontmic_0deg)
        hold on;
        plot((0:(length(rir_right_frontmic_0deg)-1))  .* (1/fs_output), rir_right_frontmic_0deg)
        plot((0:(length(rir_left_middlemic_0deg)-1))  .* (1/fs_output), rir_left_middlemic_0deg)
        plot((0:(length(rir_right_middlemic_0deg)-1)) .* (1/fs_output), rir_right_middlemic_0deg)
        plot((0:(length(rir_left_rearmic_0deg)-1))    .* (1/fs_output), rir_left_rearmic_0deg)
        plot((0:(length(rir_right_rearmic_0deg)-1))   .* (1/fs_output), rir_right_rearmic_0deg)
        xlabel('Time [sec]')
        legend('Left Front', 'Right Front', 'Left Middle', 'Right Middle', 'Left Front', 'Right Front')
        title(sprintf('Real Measured RIRs (HRIR %s)', rir_memos{rir_num}), 'Interpreter', 'none')
        
        subplot(1,2,2)
        plot((0:(length(edc_1)-1)) .* (1/fs_output), 10*log10(edc_1));
        hold on;
        plot((0:(length(edc_2)-1)) .* (1/fs_output), 10*log10(edc_2));
        plot((0:(length(edc_3)-1)) .* (1/fs_output), 10*log10(edc_3));
        plot((0:(length(edc_4)-1)) .* (1/fs_output), 10*log10(edc_4));
        plot((0:(length(edc_5)-1)) .* (1/fs_output), 10*log10(edc_5));
        plot((0:(length(edc_6)-1)) .* (1/fs_output), 10*log10(edc_6));
        grid on;
        ylim([-65 6])
        xlabel('Time [sec]')
        ylabel('dB')
        legend('Left Front', 'Right Front', 'Left Middle', 'Right Middle', 'Left Front', 'Right Front')
        title(sprintf('Energy Decay Curve (HRIR %s)', rir_memos{rir_num}), 'Interpreter', 'none')

        saveas(gcf, sprintf('RIR_EDC_%s.fig', rir_memos{rir_num}));
    end

    rir_num_4mic = rir_num_4mic + 1;
end

%% HRIR Noise


noise_memos       = {'office_ventilation',  'office_typing', 'office_telephone', 'courtyard_env', 'cafeteria_babble', 'cafeteria_aggressive_babble'};
path_preambles    = {'ambient_sound_office_II/ventilation/office_II_ventilation_1_5-00_min_bte',     ...
                     'ambient_sound_office_II/typing/office_II_typing_k1_1_3-00_min_bte',            ...
                     'ambient_sound_office_II/telephone/office_II_telephone_ringing_1_1-00_min_bte', ...
                     'ambient_sound_courtyard/courtyard_1_23-31_min_bte',                            ...
                     'ambient_sound_cafeteria/cafeteria_1_10-37_min_bte',                            ...
                     'ambient_sound_cafeteria/cafeteria_babble_1_3-23_min_bte'};

for noise_num = 1:length(noise_memos)
    [S_noise_LR_front,  fs_noise] = audioread(sprintf('HRIR_Universitat_Oldenburg/%s_front.wav', path_preambles{noise_num}));
    [S_noise_LR_middle,        ~] = audioread(sprintf('HRIR_Universitat_Oldenburg/%s_middle.wav', path_preambles{noise_num}));
    [S_noise_LR_rear,          ~] = audioread(sprintf('HRIR_Universitat_Oldenburg/%s_rear.wav', path_preambles{noise_num}));
    
    % match lengths
    min_length = min([length(S_noise_LR_front) length(S_noise_LR_middle) length(S_noise_LR_rear)]);
    S_noise_LR_front  = S_noise_LR_front(1:min_length, :);
    S_noise_LR_middle = S_noise_LR_middle(1:min_length, :);
    S_noise_LR_rear   = S_noise_LR_rear(1:min_length, :);
    
    % set desired sample rate
    if fs_noise ~= fs_output
        S_noise_LR_front  = resample(S_noise_LR_front,  fs_output, fs_noise);
        S_noise_LR_middle = resample(S_noise_LR_middle, fs_output, fs_noise);
        S_noise_LR_rear   = resample(S_noise_LR_rear,   fs_output, fs_noise);
    end

    % Collect data
    rir_databases.hrir.noise_data.noise_list{noise_num}  = [S_noise_LR_front S_noise_LR_middle S_noise_LR_rear];
    rir_databases.hrir.noise_data.noise_memos{noise_num} = noise_memos{noise_num};

    % Save 4-mic subset to a combined structures
    rir_databases_4micCombined.noise_data.noise_list{noise_num_4mic}  = [S_noise_LR_front S_noise_LR_rear];
    rir_databases_4micCombined.noise_data.noise_memos{noise_num_4mic} = noise_memos{noise_num};

    noise_num_4mic = noise_num_4mic + 1;
end


%% MYRiAD Database

rooms              = {'SAL'  };
speaker_locs_0deg  = {'S0_1' };
speaker_locs_90deg = {'S90_1'};
rir_memos          = {"SAL"  };
T60_msec_list      = {2100   };
fs_MYRiAD = 44100;

for rir_num = 1:length(rir_memos)

    speaker_loc_0deg = speaker_locs_0deg(rir_num);
    speaker_loc_90deg = speaker_locs_90deg(rir_num);

    MYRiAD_0deg  = my_load_audio_data(rooms{rir_num}, speaker_loc_0deg);
    MYRiAD_90deg = my_load_audio_data(rooms{rir_num}, speaker_loc_90deg);
    MYRiAD_audio_0deg  = MYRiAD_0deg.data{1};
    MYRiAD_audio_90deg = MYRiAD_90deg.data{1};
    
    if fs_MYRiAD ~= fs_output
        [resample_p, resample_q] = rat(fs_output / fs_MYRiAD);
        MYRiAD_audio_0deg  = resample(MYRiAD_audio_0deg, resample_p, resample_q);
        MYRiAD_audio_90deg = resample(MYRiAD_audio_90deg, resample_p, resample_q);
    end

    % Time-align RIRs and
    % Remove leading measurement noise manually (dont want to convolve this, it isnt really part of the RIR)
    [~, M] = size(MYRiAD_audio_0deg);
    tofs_est_0deg  = zeros(M,1);
    tofs_est_90deg = zeros(M,1);
    for ch = 1:M
        [MYRiAD_audio_0deg(:,ch), tofs_est_0deg(ch)]   = remove_air_meas_noise(MYRiAD_audio_0deg(:,ch));
        [MYRiAD_audio_90deg(:,ch), tofs_est_90deg(ch)] = remove_air_meas_noise(MYRiAD_audio_90deg(:,ch));
    
        % Compensate delay --> TODO Should only happen in the DAP algo i think (dont wanna get rid of ITDs)
        % tdf = [zeros(tofs_est(ch), 1) ; 1];
        % MYRiAD_audio(:,ch) = filter(tdf, 1, MYRiAD_audio(:,ch));
    end
    
    rir_left_frontmic_0deg  = MYRiAD_audio_0deg(:,1);
    rir_right_frontmic_0deg = MYRiAD_audio_0deg(:,2);
    rir_left_rearmic_0deg   = MYRiAD_audio_0deg(:,3);
    rir_right_rearmic_0deg  = MYRiAD_audio_0deg(:,4);

    rir_left_frontmic_90deg  = MYRiAD_audio_90deg(:,1);
    rir_right_frontmic_90deg = MYRiAD_audio_90deg(:,2);
    rir_left_rearmic_90deg   = MYRiAD_audio_90deg(:,3);
    rir_right_rearmic_90deg  = MYRiAD_audio_90deg(:,4);
    
    rir_databases.myriad.rir_data.rir_list_0deg{rir_num} = [rir_left_frontmic_0deg    ...
                                                       rir_right_frontmic_0deg  ...
                                                       rir_left_rearmic_0deg    ...
                                                       rir_right_rearmic_0deg ];

    rir_databases.myriad.rir_data.rir_list_90deg{rir_num} = [rir_left_frontmic_90deg    ...
                                                       rir_right_frontmic_90deg  ...
                                                       rir_left_rearmic_90deg    ...
                                                       rir_right_rearmic_90deg ];

    rir_databases.myriad.rir_data.rir_memos(rir_num)      = rir_memos{rir_num};
    rir_databases.myriad.rir_data.T60_msec_list(rir_num)  = T60_msec_list{rir_num};
    rir_databases.myriad.rir_data.tof_est_0deg{rir_num}   = tofs_est_0deg;
    rir_databases.myriad.rir_data.tof_est_90deg{rir_num}  = tofs_est_90deg;

    % Save 4-mic subset to a combined structures
    rir_databases_4micCombined.rir_data.rir_list_0deg{rir_num_4mic} = [rir_left_frontmic_0deg    ...
                                                                       rir_right_frontmic_0deg  ...
                                                                       rir_left_rearmic_0deg    ...
                                                                       rir_right_rearmic_0deg ];

    rir_databases_4micCombined.rir_data.rir_list_90deg{rir_num_4mic} = [rir_left_frontmic_90deg    ...
                                                                       rir_right_frontmic_90deg  ...
                                                                       rir_left_rearmic_90deg    ...
                                                                       rir_right_rearmic_90deg ];

    rir_databases_4micCombined.rir_data.rir_memos(rir_num_4mic)      = rir_memos{rir_num};
    rir_databases_4micCombined.rir_data.T60_msec_list(rir_num_4mic)  = T60_msec_list{rir_num};
    rir_databases_4micCombined.rir_data.tof_est_0deg{rir_num_4mic}   = tofs_est_0deg;
    rir_databases_4micCombined.rir_data.tof_est_90deg{rir_num_4mic}  = tofs_est_90deg;

    if enable_plots
        edc_1 = EDC(rir_left_frontmic_0deg);
        edc_2 = EDC(rir_right_frontmic_0deg);
        edc_3 = EDC(rir_left_rearmic_0deg);
        edc_4 = EDC(rir_right_rearmic_0deg);
    
        figure()
        set(gcf,'Position',fig_size)
        subplot(1,2,1)
        plot((0:(length(rir_left_frontmic_0deg)-1)) .* (1/fs_output), rir_left_frontmic_0deg)
        hold on;
        plot((0:(length(rir_right_frontmic_0deg)-1))  .* (1/fs_output), rir_right_frontmic_0deg)
        plot((0:(length(rir_left_rearmic_0deg)-1))    .* (1/fs_output), rir_left_rearmic_0deg)
        plot((0:(length(rir_right_rearmic_0deg)-1))   .* (1/fs_output), rir_right_rearmic_0deg)
        xlabel('Time [sec]')
        legend('Left Front', 'Right Front', 'Left Front', 'Right Front')
        title(sprintf('Real Measured RIRs (MYRiAD %s)', rir_memos{rir_num}), 'Interpreter', 'none')
        
        subplot(1,2,2)
        plot((0:(length(edc_1)-1)) .* (1/fs_output), 10*log10(edc_1));
        hold on;
        plot((0:(length(edc_2)-1)) .* (1/fs_output), 10*log10(edc_2));
        plot((0:(length(edc_3)-1)) .* (1/fs_output), 10*log10(edc_3));
        plot((0:(length(edc_4)-1)) .* (1/fs_output), 10*log10(edc_4));
        grid on;
        ylim([-65 6])
        xlabel('Time [sec]')
        ylabel('dB')
        legend('Left Front', 'Right Front','Left Front', 'Right Front')
        title(sprintf('Energy Decay Curve (MYRiAD %s)', rir_memos{rir_num}), 'Interpreter', 'none')

        saveas(gcf, sprintf('RIR_EDC_%s.fig', rir_memos{rir_num}));
    end
    rir_num_4mic = rir_num_4mic + 1;
end
rir_databases.myriad.channel_memos = ["BTE Front Left" "BTE Front Right" "BTE Rear Left" "BTE Rear Right"];


%% MYRiAD Noise

[S_noise_L_front, fs_noise] = audioread("MYRiAD_V2/audio/SAL/CP/BTELF_CP1.wav");

rir_databases.hrir.noise_data.channel_memos = ["BTE Front Left" "BTE Front Right" "BTE Rear Left" "BTE Rear Right"];

noise_memos       = {'SAL_stationary_30deg',  'SAL_music_30deg', 'SAL_cocktail_1', 'SAL_cocktail_2', 'SAL_cocktail_3', 'SAL_cocktail_4', 'SAL_cocktail_5', 'SAL_cocktail_6'};
path_preambles    = {'S30_1',                 'S30_1',    'CP',             'CP',             'CP',             'CP',             'CP',             'CP'};
path_postambles   = {'SN.wav',                'PI.wav',   'CP1.wav',        'CP2.wav',        'CP3.wav',        'CP4.wav',        'CP5.wav',        'CP6.wav'};

for noise_num = 1:length(noise_memos)
    [S_noise_L_front,  fs_noise] = audioread(sprintf("MYRiAD_V2/audio/SAL/%s/BTELF_%s", path_preambles{noise_num}, path_postambles{noise_num}));
    [S_noise_L_rear,          ~] = audioread(sprintf("MYRiAD_V2/audio/SAL/%s/BTELB_%s", path_preambles{noise_num}, path_postambles{noise_num}));
    [S_noise_R_front,         ~] = audioread(sprintf("MYRiAD_V2/audio/SAL/%s/BTERF_%s", path_preambles{noise_num}, path_postambles{noise_num}));
    [S_noise_R_rear,          ~] = audioread(sprintf("MYRiAD_V2/audio/SAL/%s/BTERB_%s", path_preambles{noise_num}, path_postambles{noise_num}));
    
    % match lengths
    min_length = min([length(S_noise_L_front) length(S_noise_L_rear) length(S_noise_R_front) length(S_noise_R_rear)]);
    S_noise_L_front  = S_noise_L_front(1:min_length);
    S_noise_L_rear   = S_noise_L_rear(1:min_length);
    S_noise_R_front  = S_noise_R_front(1:min_length);
    S_noise_R_rear   = S_noise_R_rear(1:min_length);
    
    % set desired sample rate
    if fs_noise ~= fs_output
        S_noise_L_front = resample(S_noise_L_front, fs_output, fs_noise);
        S_noise_L_rear  = resample(S_noise_L_rear,  fs_output, fs_noise);
        S_noise_R_front = resample(S_noise_R_front, fs_output, fs_noise);
        S_noise_R_rear  = resample(S_noise_R_rear,  fs_output, fs_noise);
    end

    % Collect data
    rir_databases.myriad.noise_data.noise_list{noise_num}  = [S_noise_L_front S_noise_R_front S_noise_L_rear S_noise_R_rear];
    rir_databases.myriad.noise_data.noise_memos{noise_num} = noise_memos{noise_num};

    % Save 4-mic subset to a combined structures
    rir_databases_4micCombined.noise_data.noise_list{noise_num_4mic}  = [S_noise_L_front S_noise_R_front S_noise_L_rear S_noise_R_rear];
    rir_databases_4micCombined.noise_data.noise_memos{noise_num_4mic} = noise_memos{noise_num};

    noise_num_4mic = noise_num_4mic + 1;
end

%% Save Data

save("rir_databases.mat", "rir_databases", '-v7.3'); % MAT v7.3 for large matrix
save("rir_databases_4micCombined.mat", "rir_databases_4micCombined", '-v7.3'); % MAT v7.3 for large matrix

rir_databases_4micCombined_rir_only.rir_data = rir_databases_4micCombined.rir_data;
rir_databases_4micCombined_rir_only.channel_memos = rir_databases_4micCombined.channel_memos;
rir_databases_4micCombined_rir_only.fs = rir_databases_4micCombined.fs;
save("rir_databases_4micCombined_rirOnly.mat", "rir_databases_4micCombined_rir_only");

rir_databases_4micCombined_1noise = rir_databases_4micCombined_rir_only;
rir_databases_4micCombined_1noise.noise_data.noise_list{1} = rir_databases_4micCombined.noise_data.noise_list{1};
rir_databases_4micCombined_1noise.noise_data.noise_memos{1} = rir_databases_4micCombined.noise_data.noise_memos{1};
save("rir_databases_4micCombined_1noise.mat", "rir_databases_4micCombined_1noise");


rir_databases_4micCombined_2noise = rir_databases_4micCombined_1noise;
rir_databases_4micCombined_2noise.noise_data.noise_list{2} = rir_databases_4micCombined.noise_data.noise_list{5};
rir_databases_4micCombined_2noise.noise_data.noise_memos{2} = rir_databases_4micCombined.noise_data.noise_memos{5};
save("rir_databases_4micCombined_2noise.mat", "rir_databases_4micCombined_2noise");
