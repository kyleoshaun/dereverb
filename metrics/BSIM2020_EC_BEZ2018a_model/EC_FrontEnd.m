function [mixed_output] = EC_FrontEnd(teststim_left, teststim_right, Fs_stim)

%% Add Paths

% Calling script must include:
% addpath('../BSIM_2020')
% addpath('../BSIM_2020/SRMRToolbox-master')
% addpath('../BSIM_2020/SRMRToolbox-master/libs/')
% addpath('../BSIM_2020/SRMRToolbox-master/libs/vad/')
% addpath('../BSIM_2020/SRMRToolbox-master/libs/PreProc/')
% addpath('../BSIM_2020/SRMRToolbox-master/libs/Gammatonegram/')
% addpath('../BSIM_2020/SRMRToolbox-master/libs/Auditory/')
% 
% addpath(genpath('../BSIM_2020/Anechoic/'))
% addpath('../BSIM_2020/functions/')
% addpath('../BSIM_2020/functions/gammatonefilter/')
% addpath('../BSIM_2020/functions/ltfat/')

%% Init EC (BSIM 2020)

%Model_init;

%----------------Model parameters-----------------------------------------%
% Used defaults from BSIM20()
model_param.fs = Fs_stim; % sampling rate 
model_param.fmin = 125; % lowest filter of gammatone filterbank
model_param.fmax = 8000; % highest filter of gammatone filterbank --> Bruce et al Model runs at 100 kHz, but we dont actually need 50kHz BW
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


%% Run EC Front-End

% TEMP
%teststim_left = resample(teststim_left,model_param.fs,Fs_stim);
%refstim = resample(refstim,model_param.fs,Fs_stim);
% END

if length(teststim_left) ~= length(teststim_right)
    error("Left and Right Signal lengths must match")
end

mixed_input = [teststim_left teststim_right];
%out_struct = BSIM20('RequiredSignal',mixed_input);
out_struct = BSIM20('RequiredSignal',mixed_input, 'model_params',model_param);

teststim = out_struct.signals.MixsigSyn'; %out_struct.signals.MixsigSyn;

mixed_output = teststim; % Output for reference

end