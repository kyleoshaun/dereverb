addpath('../../BEZ2018a_model/')

% Check to see if running under Matlab or Octave
if exist ('OCTAVE_VERSION', 'builtin') ~= 0
  pkg load signal;
  if exist('rms')<1
    rms = @(x) sqrt(mean(x.^2));
  end
end

if exist('parfor','builtin') % check if the Matlab Parallel Computation
                             % Toolbox is installed and use appropriate
                             % function
    generate_neurogram_function = @generate_neurogram_BEZ2018a_stmi_parallelized;
    disp('Using parallelized version of neurogram generation function')
else
    generate_neurogram_function = @generate_neurogram_BEZ2018a_stmi;
    disp('Using serial version of neurogram generation function')
end

% Set "reference" audiogram for STMI
ag_fs_ref = [0 20e3];
ag_dbloss_ref = [0 0]; % Normal hearing

% Set "test" audiogram for STMI
ag_fs_test = [0 20e3];
ag_dbloss_test = [0 0]; % Normal hearing

species = 2; % human model with tuning based on Shera et al. (2002)

rv = 2.^(1:0.5:5);  % cortical (temporal) modulation filter rates
sv = 2.^(-2:0.5:3); % cortical (spectral) modulation filter scales

[refstim, Fs_stim] = audioread('defineit.wav');

stimdb = 65; % speech level in dB SPL
SNR = 5; % in dB
% SNR = inf; % in dB; inf -> no background noise

refstim = refstim/rms(refstim)*20e-6*10^(stimdb/20);
noisestim = randn(size(refstim))*rms(refstim)*10^(-SNR/20);

teststim = refstim + noisestim;

[neurogram_stmi_ref,t_stmi,CFs] = generate_neurogram_function(refstim,Fs_stim,species,ag_fs_ref,ag_dbloss_ref);
base_sig_ref = generate_base_signal(refstim);
[neurogram_base_ref,t_stmi,CFs] = generate_neurogram_function(base_sig_ref,Fs_stim,species,ag_fs_ref,ag_dbloss_ref);

[neurogram_stmi_test,t_stmi,CFs] = generate_neurogram_function(teststim,Fs_stim,species,ag_fs_test,ag_dbloss_test);
base_sig_test = generate_base_signal(teststim);
[neurogram_base_test,t_stmi,CFs] = generate_neurogram_function(base_sig_test,Fs_stim,species,ag_fs_test,ag_dbloss_test);

response_endtime = length(refstim)/Fs_stim+10e-3; % Calculate end time of the neural response to the stimulus
[tmp, ind_stmi] = min(abs(t_stmi-response_endtime)); % Find the corresponding time index in the stmi neurogram

% Generate the cortical responses to the reference signal and reference base signal neurograms and
% compute the reference "template" T matrix
a1_stmi_ref = abs(ngram2cortex(neurogram_stmi_ref(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
a1_base_ref = abs(ngram2cortex(neurogram_base_ref(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
T  = max(a1_stmi_ref - a1_base_ref,0);
TT = sum(sum(T(:).^2));

% Generate the cortical responses to the test signal and test base signal neurograms and
% compute the test "noisy" N matrix
a1_stmi_test = abs(ngram2cortex(neurogram_stmi_test(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
a1_base_test = abs(ngram2cortex(neurogram_base_test(:,1:ind_stmi)',diff(t_stmi(1:2)),rv,sv)); % Neurogram needs to be transposed so the response for each CF is in a column -> 128 columns
N  = max(a1_stmi_test - a1_base_test,0);
NN = sum(sum((max(T(:)-N(:),0)).^2));

STMI = 1-NN/TT

ng1=figure;
set(ng1,'renderer','painters');
winlen = 256; % Window length for the spectrogram analyses

sp1 = subplot(2,2,1);
[s,f,t] = specgram([refstim; eps*ones(round(t_stmi(end)*Fs_stim)-length(refstim),1)],winlen,Fs_stim,winlen,0.25*winlen);
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
xlabel('Time');
ylabel('Frequency (kHz)');
title('Reference Stimulus Spectrogram')
xl = xlim;

sp2 = subplot(2,2,2);
[s,f,t] = specgram([teststim; eps*ones(round(t_stmi(end)*Fs_stim)-length(teststim),1)],winlen,Fs_stim,winlen,0.25*winlen);
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Test Stimulus Spectrogram')
xlim(xl)

sp3=subplot(2,2,3);
plot_neurogram(t_stmi,CFs,neurogram_stmi_ref,sp3);
caxis([0 200])
title('STMI Neurogram for Reference Stimulus')
xlim(xl)

sp4=subplot(2,2,4);
plot_neurogram(t_stmi,CFs,neurogram_stmi_test,sp4);
caxis([0 200])
title(['STMI Neurogram for Test Stimulus: STMI = ' num2str(STMI)])
xlim(xl)

if exist ('OCTAVE_VERSION', 'builtin') == 0
drawnow;
ng1.WindowState = 'Maximized';
end
