close all
clear
clc

addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/')
addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/')
addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/courtyard/')
addpath('../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/cafeteria/')
addpath('../../RIR_Databases/MYRiAD_V2/tools/MATLAB/')

%% HRIR

rooms = cell(3,1);
rooms{1} = 'Office_II';
rooms{2} = 'Courtyard';
rooms{3} = 'Cafeteria';

for room_num = 1:length(rooms)

    room             = rooms{room_num};
    head_orientation = 1;
    speaker_loc      = 'A';
    data_set         = 'bte';
    HRIR_data = loadHRIR(room, head_orientation, speaker_loc, data_set);
    fs = HRIR_data.fs;
    
    b_rir_1 = HRIR_data.data(:,1);
    b_rir_2 = HRIR_data.data(:,2);
    b_rir_3 = HRIR_data.data(:,3);
    b_rir_4 = HRIR_data.data(:,4);
    b_rir_5 = HRIR_data.data(:,5);
    b_rir_6 = HRIR_data.data(:,6);
    
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

end


%% MYRiAD

    MYRiAD_audio = my_load_audio_data;

    fs = 44100;
    
    b_rir_1 = MYRiAD_audio.data{1}(:,1);
    b_rir_2 = MYRiAD_audio.data{1}(:,2);
    b_rir_3 = MYRiAD_audio.data{1}(:,3);
    b_rir_4 = MYRiAD_audio.data{1}(:,4);
    
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