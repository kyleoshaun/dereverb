addpath('../../BEZ2018a_model/')

% Check to see if running under Matlab or Octave
if exist ('OCTAVE_VERSION', 'builtin') ~= 0
  pkg load signal;
  if exist('rms','builtin')<1
    rms = @(x) sqrt(mean(x.^2));
  end
end

if exist('parfor','builtin') % check if the Matlab Parallel Computation
                             % Toolbox is installed and use appropriate
                             % function
    generate_neurogram_function = @generate_neurogram_BEZ2018a_parallelized;
    disp('Using parallelized version of neurogram generation function')
else
    generate_neurogram_function = @generate_neurogram_BEZ2018a;
    disp('Using serial version of neurogram generation function')
end

% Set "reference" audiogram
ag_fs_ref = [0 20e3];
ag_dbloss_ref = [0 0]; % Normal hearing

% Set "test" audiogram
ag_fs_test = [125 250 500 1e3 2e3 4e3 8e3];
ag_dbloss_test = [0 0 0 0 0 0 0]; % Normal hearing
%ag_dbloss_test = [0 0 0 20 40 60 80]; % Example high-freq hearing loss



species = 2;

% NSIM parameters
weights.alpha = 1.0;
weights.beta = 0.0;
weights.gamma = 1.0;
window_type = 0;
win3x3 = ones(3,3);

[refstim, Fs_stim] = audioread('defineit.wav');

stimdb = 65; % speech level in dB SPL
SNR = 0; % in dB
% SNR = inf; % in dB; inf -> no background noise

refstim = refstim/rms(refstim)*20e-6*10^(stimdb/20);
noisestim = randn(size(refstim))*rms(refstim)*10^(-SNR/20);

teststim = refstim + noisestim;

[neurogram_ft_ref,neurogram_mr_ref,neurogram_Sout_ref,t_ft,t_mr,t_Sout,CFs] = generate_neurogram_function(refstim,Fs_stim,species,ag_fs_ref,ag_dbloss_ref);

[neurogram_ft_test,neurogram_mr_test,neurogram_Sout_test,t_ft,t_mr,t_Sout,CFs] = generate_neurogram_function(teststim,Fs_stim,species,ag_fs_test,ag_dbloss_test);

response_endtime = length(refstim)/Fs_stim+10e-3; % Calculate end time of the neural response to the stimulus
[tmp, ind_mr] = min(abs(t_mr-response_endtime)); % Find the corresponding time index in the mr neurogram
[tmp, ind_ft] = min(abs(t_ft-response_endtime)); % Find the corresponding time index in the ft neurogram

% Scaling method used by Hines and Harte (Speech Comm 2010, 2012)
% scl_mr = 255/max(max(neurogram_mr_ref(:,1:ind_mr)));
% scl_ft = 255/max(max(neurogram_ft_ref(:,1:ind_ft)));

% New scaling method developed by M. R. Wirtzfeld (see Wirtzfeld et al., JARO 2017)
scl_mr = 1/50/t_mr(2);
scl_ft = 1/50/t_ft(2);

[NSIM_MR, ssim_mr] = mssim_v2(scl_mr*neurogram_mr_ref(:,1:ind_mr), scl_mr*neurogram_mr_test(:,1:ind_mr), weights, win3x3, window_type );

[NSIM_FT, ssim_ft] = mssim_v2(scl_ft*neurogram_ft_ref(:,1:ind_ft), scl_ft*neurogram_ft_test(:,1:ind_ft), weights, win3x3, window_type );

NSIM_MR

NSIM_FT

%%

ng1=figure;
set(ng1,'renderer','painters');
winlen = 256; % Window length for the spectrogram analyses

sp1 = subplot(3,2,1);
[s,f,t] = specgram([refstim; eps*ones(round(t_mr(end)*Fs_stim)-length(refstim),1)],winlen,Fs_stim,winlen,0.25*winlen);
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

sp2 = subplot(3,2,2);
[s,f,t] = specgram([teststim; eps*ones(round(t_mr(end)*Fs_stim)-length(teststim),1)],winlen,Fs_stim,winlen,0.25*winlen);
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

sp3=subplot(3,2,3);
plot_neurogram(t_mr,CFs,neurogram_mr_ref,sp3);
caxis([0 80])
title('Mean-rate Neurogram for Reference Stimulus')
xlim(xl)

sp4=subplot(3,2,4);
plot_neurogram(t_mr,CFs,neurogram_mr_test,sp4);
caxis([0 80])
title(['Mean-rate Neurogram for Test Stimulus: NSIM\_MR = ' num2str(NSIM_MR)])
xlim(xl)

sp5=subplot(3,2,5);
plot_neurogram(t_ft,CFs,neurogram_ft_ref,sp5);
caxis([0 20])
title('Fine-timing Neurogram for Reference Stimulus')
xlim(xl)

sp6=subplot(3,2,6);
plot_neurogram(t_ft,CFs,neurogram_ft_ref,sp6);
caxis([0 20])
title(['Fine-timing Neurogram for Test Stimulus: NSIM\_FT = ' num2str(NSIM_FT)])
xlim(xl)

