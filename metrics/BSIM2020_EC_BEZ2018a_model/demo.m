close all
clear
clc

%% NEW REVERB DEMO

addpath('../../dereverb/LPC/../../samples')
addpath('../../dereverb/LPC/../../RIR_Databases/AIR_1_4_BinauralDatabase/')
addpath('../../dereverb/LPC/../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/')
addpath('../../dereverb/LPC/../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/office_II/')
addpath('../../dereverb/LPC/../../RIR_Databases/HRIR_Universitat_Oldenburg/HRIR_database_mat/hrir/courtyard/')

[refstim, Fs_stim] = audioread("SA1.wav");

% TODO: Scale to 65 dB SPL AFTER BRIR application (SPL at ear drum not source)
% -- Also take a look at effect of EC

stimdb = 65; % speech level in dB SPL 
% SNR = 0; % in dB
% % SNR = inf; % in dB; inf -> no background noise

refstim = refstim/rms(refstim)*20e-6*10^(stimdb/20);
% noisestim = randn(size(refstim))*rms(refstim)*10^(-SNR/20);

% Option 2: Real AIR
% Channels:
% - 1: Left Front
% - 2: Right Front
% - 3: Left Middle
% - 4: Right Middle
% - 5: Left Rear
% - 6: Front Rear
head_orientation = 2;
speaker_loc      = 'A';
data_set         = 'bte';
HRIR_data = loadHRIR('office_II', head_orientation, speaker_loc, data_set);
if HRIR_data.fs ~= Fs_stim
    [resample_p, resample_q] = rat(Fs_stim / HRIR_data.fs);
    HRIR_data.data = resample(HRIR_data.data, resample_p, resample_q);
end

b_air_1 = HRIR_data.data(:,1); % LEFT FRONT
b_air_2 = HRIR_data.data(:,2); % RIGHT FRONT

teststim_left  = filter(b_air_1, 1, refstim);
teststim_right = filter(b_air_2, 1, refstim);

%% OLD NOISE DEMO

% %Fs_stim = 100000;
% [refstim, Fs_stim] = audioread('defineit.wav');
% 
% stimdb = 65; % speech level in dB SPL
% SNR = 0; % in dB
% % SNR = inf; % in dB; inf -> no background noise
% 
% refstim = refstim/rms(refstim)*20e-6*10^(stimdb/20);
% noisestim = randn(size(refstim))*rms(refstim)*10^(-SNR/20);
% 
% teststim = refstim + noisestim;
% 
% teststim_left  = teststim;
% teststim_right = teststim;

%% RUN


[STMI, NSIM_FT, ec_out] = binaural_nsim_and_stmi(teststim_left, teststim_right, refstim, Fs_stim, stimdb);

fprintf("RESULTS: --------------------------------------\n")
fprintf("STMI    = %.8f\n", STMI);
fprintf("NSIM_FT = %.8f\n", NSIM_FT);
fprintf("-----------------------------------------------")