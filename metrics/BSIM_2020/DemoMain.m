% Usage: Type <DemoMain> in the command line
% Framework to test the binaural speech intelligibility model BSIM2020 
% published in "Trends in Hearing - complete journal stuff here!!!CH"
% to predict the SRTs from the Anechoic Experiment of Beutelmann and Brand 
% (2006). In this experiment, SRTs are simulated for speech located at 0° in
% the horizontal plane, and noise located at different positions [angles].
% The binaural processing in this model works blindly and only requires the
% mixture of speech an noise [required signal]. The back-end used here is the SII. 
% SII values are produced for different signal-to-noise ratios. Depending
% on the back-end you want to use, optional signals have to be used. They
% are processed in the same way as the required signal. In this Demo, the
% optional signals are clean speech and noise. 
% This release uses the gammatone filters from the AMToolbox.
%
% Authors: Christopher Hauth <christopher.hauth@uni-oldenburg.de>
%          Dr. Thomas Brand  <thomas.brand@uni-oldenburg.de>
% Date : 22.10.2020
% Version: 1.0
%-------------------------------------------------------------------------%
clear;
close all;
clc;
% Results can be saved in the Results Folder
if ~exist('Results','Dir')
    mkdir Results;
end

if exist('ini_flag','var')
    disp('Initialization done!');
else
    Model_init;
    ini_flag = 1;
end

list_clean = dir('Anechoic/0');
%---------------Experimental Conditions-----------------------------------%
% Define your experimental conditions, number of Monte Carlo simulations
% (binaural processing) and number of sentences (statistics across sentences)
% For Matrix type sentences it is recommended to use 10 sentences, where each word
% of the test appears once.
iNumofSentences = 1;
iNumofMonteCarlo = 10;

sentences_clean = {list_clean.name};
sentence_choose_clean  = sentences_clean(3:3+iNumofSentences);%choose your sentences (needs to be precise for publication, when we know the number of sentences available [Comment: CH])
randomizer = randperm(length(sentence_choose_clean));

vangles_test  = [-140 -100 -45 0 45 80 125 180] ; % angles from which the noise is presented
% Vector of input SNRs. For each SNR, (iNumofSentences x iNumofMonteCarlo) SII values are obtained.
% Make sure to test different different SNRs in order to be able to map the
% SII to an SRT, e.g.: 
 vSNR_test = -20:0;

% If you want to know the SII of a single SNR, use only one value:
 %vSNR_test = -18;                                
%----------------Model parameters-----------------------------------------%
model_param.fs = 44100; % sampling rate
model_param.fmin = 150; % lowest filter of gammatone filterbank
model_param.fmax = 8500; % highest filter of gammatone filterbank
model_param.f_target = 500; % specified frequeny that will have a matched filter
model_param.bin_err = 1; % enable(1)/disable(0) binaural processing inaccuracies (important for Monte-Carlo Simulations)
model_param.ERB_factor = 1;% define bandwidth of filters
model_param.long_term = 1; % use long-term model (1) or short-term model (0)
model_param.Filterbank = 'GT';% define filtering: GT is gammatone filter
% Model parameters don't have to be specified. These are also the default values,
% which are used in the model
% Please note that the modified SRMR model, which performs the channel
% selection either considers the whole signal as a single time frame (long
% term or uses the same time constants as the EC processing and BE
% processing if the short time model is used.
%-------------------------------------------------------------------------%
%% Calibration
% the calibration factor is mean level between the ears
% The calibration can be adjusted for your needs. For the SII, which is used here, 65
% dB FS (relative to full scale) is assumed to be 65 dB SPL.
% However, if you only aim for the
% resynthesized output, please adjust this level to avoid clipping. 
lev_desired = 65;   % the value 65 is required here for correct use of the SII;

% Use the co-located noise condition to calibrate the input (For the speeech signal, also the noise is used)
[calibnoise fs_n]= audioread(sprintf('Anechoic/0/%d_speaker_reference_olnoise.wav',0));
calibnoise = calibnoise(1:end-round(1.5*fs_n),:);
    
lev_Speech = 20*log10(rms(calibnoise));  % actual rms-level of speech
lev_S = mean(lev_Speech);                % the mean between both ears is considered
Delta_L_speech = lev_desired - lev_S;    % Calibration Gain is the difference between the desired level and the actual level
Delta_L_speech_lin = 10.^((Delta_L_speech)/20); % Convert to linear gain

% Similar calibration of the noise
lev_Noise = 20*log10(rms(calibnoise));  % actual rms-level of the noise
lev_Noise = mean(lev_Noise);                % reference is MEAN level between the two ears
Delta_L_noise = lev_desired - lev_Noise;    % Calibration Gain is the difference between the desired level and the actual level
Delta_L_noise_lin = 10.^((Delta_L_noise)/20);% Convert to linear gain
%-------------------------------------------------------------------------%
% Iterate through all angles of the noise position
for mm = 1:length(vangles_test)
    % Iterate through all SNRs
    for kk = 1:length(vSNR_test)
        % Iterate through the different sentences
        for ll = 1:iNumofSentences
            sentence_clean = sentence_choose_clean{randomizer(ll)};         
            % Read sentences and noise from wav files
            [speech_clean fs_s]= audioread(['Anechoic/0/' sentence_clean]);
            [noise fs_n]= audioread(sprintf('Anechoic/%d/%d_speaker_reference_olnoise.wav',vangles_test(mm),vangles_test(mm)));
            
            % resample signals if necessary 
            speech_clean = resample(speech_clean,model_param.fs,fs_s);
            noise = resample(noise,model_param.fs,fs_n);
 
            % Get length of the speech signal
            lenSpeech = length(speech_clean);
            lenNoise = length(noise);
            
            % Truncate noise to have the same lenght as speech
            noise = noise(1:lenSpeech,:);
            speech_clean = speech_clean(1:lenSpeech,:);
            % adjust level of speech:
            speech_clean = Delta_L_speech_lin.*speech_clean;
            
            % adjust level of noise:
            noise = Delta_L_noise_lin.*noise;
            
            % adjust SNR of mixed input signal (speech + noise)
            % This is a required signal:
            mixed_input = 10.^((vSNR_test(kk))/20).*speech_clean+noise;
            inputLen = length(mixed_input);
            
            % Adjust level of the clean speech if you want to 
            % use it as an optional input:
            speech_clean_proc = 10.^((vSNR_test(kk))/20).*speech_clean;
            
            % All optional signals are arranged in a matrix.
            % Here: [S_l(:) S_r(:) N_l(:) N_r(:)]
            OptionalSignals = [speech_clean_proc noise];
            
            % Apply the binaural model to the mixed signal 
            % (and optionally to the clean speech and noise)
            % Monte Carlo simulations are used to model the binaural uncerntainty
            for oo=1:iNumofMonteCarlo
                % Do binaural processing:
                % out_struct contains the processed mixed signal as well as
                % the processed optional signals: Moreover, as the SII is
                % used as back-end, it also contains the frequency-specific 
                % levels

                 out_struct = BSIM20('RequiredSignal',mixed_input,...
                  'OptionalSignal',OptionalSignals,'model_params',model_param);

                % Example of BSIM2020 without any optional signals. Here, only
                % the processed mixed signals are in the output struct 
                 out_struct2 = BSIM20('RequiredSignal',mixed_input,...
                 'model_params',model_param);

                % Speech Intelligibility back-end: (SII in this example)
                % Use your speech intelligibility back-end here
                [sii_min_temp(oo),A,Z] = SII(out_struct.levels.LevelOptSig1min,out_struct.levels.LevelOptSig2min,-Inf*ones(30,1),2,0); 
                [sii_max_temp(oo),A,Z] = SII(out_struct.levels.LevelOptSig1max,out_struct.levels.LevelOptSig2max,-Inf*ones(30,1),2,0); 
                [sii_syn_temp(oo),A,Z] = SII(out_struct.levels.LevelOptSig1syn,out_struct.levels.LevelOptSig2syn,-Inf*ones(30,1),2,0); 
                [sii_L_temp(oo),A,Z]   = SII(out_struct.levels.LevelOptSig1L,out_struct.levels.LevelOptSig2L,-Inf*ones(30,1),2,0); 
                [sii_R_temp(oo),A,Z]   = SII(out_struct.levels.LevelOptSig1R,out_struct.levels.LevelOptSig2R,-Inf*ones(30,1),2,0); 
                          
                sii_min_all(ll,kk,oo) = sii_min_temp(oo);
                sii_max_all(ll,kk,oo) = sii_max_temp(oo);
                sii_syn_all(ll,kk,oo) = sii_syn_temp(oo);
                sii_L_all(ll,kk,oo)   = sii_L_temp(oo);
                sii_R_all(ll,kk,oo)   = sii_R_temp(oo);     
            end
        end
    end
    % Save the resulting SII-values:
     save(sprintf('Results/SII_All_Anechoic_%d_Min.mat',vangles_test(mm)),'sii_min_all');
     save(sprintf('Results/SII_All_Anechoic_%d_Max.mat',vangles_test(mm)),'sii_max_all');
     save(sprintf('Results/SII_All_Anechoic_%d_Syn.mat',vangles_test(mm)),'sii_syn_all');
     save(sprintf('Results/SII_All_Anechoic_%d_L.mat',vangles_test(mm)),'sii_L_all');
     save(sprintf('Results/SII_All_Anechoic_%d_R.mat',vangles_test(mm)),'sii_R_all'); 
end
%--------------------Licence ---------------------------------------------
% Copyright (c) <2020> Christopher F. Hauth
% Dept. Medical Physics and Acoustics
% Carl von Ossietzky University Oldenburg 
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to 
% permit persons to whom the Software is furnished to do so, subject 
% to the following conditions:
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% END OF FILE