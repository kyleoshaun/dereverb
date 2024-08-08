close all
clear % clear all variables from workspace
clc % clear command window

addpath('../LIVELabSpeakerLoc/')
addpath('../../../Noise_Databases/ARTE_Ambisonics/ARTE database structure/baseFunctions/')
addpath('../../../Noise_Databases/ARTE_Ambisonics/ARTE database structure/MOA31ch/noiseFiles/')

%% LIVELab

load('sphcord_mat.mat')

az_all = sphcordmat(:,2);
el_all = sphcordmat(:,3);

%idx = [13 14 23 24 33 34 47 50]';
idx_el1 = [13 14 23 24 33 34 47 50]';
idx_el2 = [];%4 5 6 7 8 9]';% 10 11 12]';
idx     = [idx_el1 ; idx_el2];

az  = az_all(idx);
el1 = el_all(idx_el1);
el2 = el_all(idx_el2);

el1 = mean(el1) * ones(length(el1),1);
el2 = mean(el2) * ones(length(el2),1);
el  = [el1 ; el2];

LsPos = [az el];


%LsPos = [(0:22.5:337.5)' zeros(16,1)]; % horizontal ring of 16 equidistant loudspeakers

LsPos = LsPos/180*pi;% transform loudspeaker locations from degrees to radians

Nlsp = size(LsPos,1);% number of loudspeakers (here 41)

plot_LIVELab(idx)

%%% -----------------------------------------------------------------------
%%% Initialization
%%% -----------------------------------------------------------------------


%addpath('baseFunctions') %set path to main MOA functions

M2 = 3;%2D HOA order - must be <= 7 (Note: requires horizontal ring of k >= 2*M2+1 loudspeakers)
M3 = 1;%3D HOA order - must be <= 4 (Note: requires regular array of k >= (M3+1)^2 loudspeakers)

Nlsp_min = (M3+1)^2 + (M2-M3)*2;
if (Nlsp < Nlsp_min)
    error('Not enough loudspeakers (%d < %d)', Nlsp, Nlsp_min)
end

%[LsPos,Nlsp] = Lsp16ch;% Creates Nlsp x 2 matrix that defines the applied loudspeaker array layout. 
                       % Here, a regular horizontal 16 channel 2D array is considered as an example.

%size_LsPos = size(LsPos);
%N_ls = size_LsPos(1);

[D,MOAch] = decMatrix(LsPos, M3, M2);% derive basic decoding matrix D (and utilized MOA channels)

Nbits = 32;% number of Bits per sample for writing multi-channel loudspeaker wav-file
fs = 44100;% Sampling frequency (Hz) NOte that all provided sound files were sampled at 44100 Hz

%dBFS = -12;
%scale = 10 ^ (dBFS/20);
gain = 4;

%%% -----------------------------------------------------------------------
%%% Decode MOA
%%% -----------------------------------------------------------------------

%N = size(sTreverbLsp,1);% length of reverberant speech (samples)

%%% Read MOA-coded noise file that correpsonds to the above RIR and
%%% has the same length as the reverberant speech.
%nameNoise = '12_Food_Court_1_MOA_31ch.wav';
nameNoise = '06_Diffuse_noise_MOA_31ch.wav';
sNoiseMOA = readMOAfile(nameNoise,MOAch);%,N);

% Decode MOA-coded noise signal into multi-channel loudspeaker signal
sNoiseLsp = MOA2Lsp(sNoiseMOA,D);

% Saving wav-file to disk
for channel_idx = 1:Nlsp
    %nameSave = sprintf('decoded/12_Food_Court_1_MOA_31ch_%d.wav', channel_idx);
    nameSave = sprintf('decoded/s_%d.wav', idx(channel_idx));
    audiowrite(nameSave,gain .* sNoiseLsp(:,channel_idx),fs,'BitsPerSample',Nbits)
end

fprintf("\nELEVATION 1 (%.2f m):\n", el1(1))
for channel_idx = 1:length(idx_el1)
    fprintf("rms(%d) = %.4f dB\n", channel_idx, db(rms(sNoiseLsp(:,channel_idx))))
end

if (isempty(el2) == false)
    fprintf("\nELEVATION 2 (%.2f m):\n", el2(1))
    for channel_idx = (length(idx_el1)+1):length(idx)
        fprintf("rms(%d) = %.4f dB\n", channel_idx, db(rms(sNoiseLsp(:,channel_idx))))
    end
end


[s_13,fs] = audioread('./decoded/s_13.wav');
[s_14,fs] = audioread('./decoded/s_14.wav');
[s_23,fs] = audioread('./decoded/s_23.wav');
[s_24,fs] = audioread('./decoded/s_24.wav');

figure()
subplot(4,1,1)
plot(s_13)
ylim([-0.02, 0.02])
subplot(4,1,2)
plot(s_14)
ylim([-0.02, 0.02])
subplot(4,1,3)
plot(s_23)
ylim([-0.02, 0.02])
subplot(4,1,4)
plot(s_24)
ylim([-0.02, 0.02])