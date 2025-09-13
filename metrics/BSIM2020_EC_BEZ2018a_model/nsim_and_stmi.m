function [STMI, NSIM_FT, NSIM_MR] = nsim_and_stmi(teststim, refstim, Fs_stim, HL, stimdb, results_dir)
% Computes Spectro-Temporal Modulation Index (STMI) and 
% Mean-Rate/Fine-Timing Neurogram Similarity Index Measure (MR NSIM / FT NSIM)
%
% Syntax:
%   [STMI, NSIM_FT, NSIM_MR] = nsim_and_stmi(teststim, refstim, Fs_stim, HL)
% 
% Inputs: 
% - teststim: Test signal in Pa at the ear drum (Distorted, processed)
% - refstim: Reference Signal in Pa (Clean, undistorted, unprocessed)
% - Fs_stim: Sample rate for stimuli (Hz)
% - HL: Hearing Loss Vector following HASPI convention
% 	    (1,6) vector of hearing loss at the 6 audiometric frequencies
%	    [250, 500, 1000, 2000, 4000, 6000] Hz.
% - stimdb: Stimulus level in dB SPL
% - results_dir: Directory for saving results
%
% Outputs:
% - STMI: Spectro-Temporal Modulation Index
% - NSIM_FT: Fine-Timing Neurogram Similarity Index
% - NSIM_MR: Mean-Rate Neurogram Similarity Index

enable_plots = true;

%% Init Model (BEZ2018a)

% Set "reference" audiogram
ag_fs_ref = [0 20e3];
ag_dbloss_ref = [0 0]; % Normal hearing

% Set "test" audiogram
ag_fs_test     = [250 500 1000 2000 4000 6000];
ag_dbloss_test = HL;

species = 2; % Human

%% Other Params

% NSIM Params
% weights.alpha = 1.0;
% weights.beta = 0.0;
% weights.gamma = 1.0;

% IDEAL WEIGHTS -- PRODUCES COMPLEX NSIM VALUES
% weights_NSIM_FT.alpha = 1.0;
% weights_NSIM_FT.beta = 0.0;
% weights_NSIM_FT.gamma = 0.65;
% weights_NSIM_MR.alpha = 0.95;
% weights_NSIM_MR.beta = 0.0;
% weights_NSIM_MR.gamma = 0.6;

% NSIM Weights
weights_NSIM_FT.alpha = 1.0;
weights_NSIM_FT.beta = 0.0;
weights_NSIM_FT.gamma = 1.0;
weights_NSIM_MR.alpha = 1.0;
weights_NSIM_MR.beta = 0.0;
weights_NSIM_MR.gamma = 1.0;


window_type = 0;
win3x3 = ones(3,3);

% STMI Params
rv = 2.^(1:0.5:5);  % cortical (temporal) modulation filter rates
sv = 2.^(-2:0.5:3); % cortical (spectral) modulation filter scales

% Parameters for STMI/NSIM Synthesis
% Regression of Fine-Timing NSIM + STMI performed by Rationalized Arcsine Transformation (RAU)
% SI_SYN = rau_b0 + rau_b1*STMI + rau_b2*NSIM_FT + rau_b3*STMI*NSIM_FT;
% rau_b0 = 0;
% rau_b1 = 1.357;
% rau_b2 = 0.444;
% rau_b3 = 0.00598;


%% NSIM

% SETUP

if exist('parfor','builtin') % check if the Matlab Parallel Computation
                             % Toolbox is installed and use appropriate
                             % function
    generate_neurogram_function = @generate_neurogram_BEZ2018a_parallelized;
    disp('Using parallelized version of neurogram generation function')
else
    generate_neurogram_function = @generate_neurogram_BEZ2018a;
    disp('Using serial version of neurogram generation function')
end

%refstim = 0;
%teststim = 0;

% NSIM CALCULATION

[neurogram_ft_ref,neurogram_mr_ref,neurogram_Sout_ref,t_ft,t_mr,t_Sout,CFs] = generate_neurogram_function(refstim,Fs_stim,species,ag_fs_ref,ag_dbloss_ref);

[neurogram_ft_test,neurogram_mr_test,neurogram_Sout_test,t_ft,t_mr,t_Sout,CFs] = generate_neurogram_function(teststim,Fs_stim,species,ag_fs_test,ag_dbloss_test);

response_endtime = length(refstim)/Fs_stim+10e-3; % Calculate end time of the neural response to the stimulus
[tmp, ind_mr] = min(abs(t_mr-response_endtime)); % Find the corresponding time index in the mr neurogram
[tmp, ind_ft] = min(abs(t_ft-response_endtime)); % Find the corresponding time index in the ft neurogram

% Scaling method used by Hines and Harte (Speech Comm 2010, 2012)
%scl_mr = 255/max(max(neurogram_mr_ref(:,1:ind_mr)));
%scl_ft = 255/max(max(neurogram_ft_ref(:,1:ind_ft)));

% New scaling method developed by M. R. Wirtzfeld (see Wirtzfeld et al., JARO 2017)
scl_mr = 1/50/t_mr(2);
scl_ft = 1/50/t_ft(2);

[NSIM_MR, ssim_mr] = mssim_v2(scl_mr*neurogram_mr_ref(:,1:ind_mr), scl_mr*neurogram_mr_test(:,1:ind_mr), weights_NSIM_MR, win3x3, window_type );

[NSIM_FT, ssim_ft] = mssim_v2(scl_ft*neurogram_ft_ref(:,1:ind_ft), scl_ft*neurogram_ft_test(:,1:ind_ft), weights_NSIM_FT, win3x3, window_type );

% NSIM PLOTS

if enable_plots
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

    if any(HL)
        sgtitle("NSIM Neurograms (HI Listener)")
        saveas(gcf, sprintf('%s/NSIM_Neurograms_HI.fig', results_dir));
    else
        sgtitle("NSIM Neurograms (NH Listener)")
        saveas(gcf, sprintf('%s/NSIM_Neurograms_NH.fig', results_dir));
    end
end


%% STMI

% STMI CALC SETUP

if exist('parfor','builtin') % check if the Matlab Parallel Computation
                             % Toolbox is installed and use appropriate
                             % function
    generate_neurogram_function = @generate_neurogram_BEZ2018a_stmi_parallelized;
    disp('Using parallelized version of neurogram generation function')
else
    generate_neurogram_function = @generate_neurogram_BEZ2018a_stmi;
    disp('Using serial version of neurogram generation function')
end

rv = 2.^(1:0.5:5);  % cortical (temporal) modulation filter rates
sv = 2.^(-2:0.5:3); % cortical (spectral) modulation filter scales

% refstim = 0;
% teststim = 0;

% STMI CALC

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

STMI = 1-NN/TT;

% STMI PLOTS

if enable_plots
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

    if any(HL)
        sgtitle("STMI Neurograms (HI Listener)")
        saveas(gcf, sprintf('%s/STMI_Neurograms_HI.fig', results_dir));
    else
        sgtitle("STMI Neurograms (NH Listener)")
        saveas(gcf, sprintf('%s/STMI_Neurograms_NH.fig', results_dir));
    end
end

%% Half-wave Rectify NSIM to account for potential negative numbers 
%  Note: Negative numbers occur somtimes in structure term of SSIM in 
%  cross-covariance term

NSIM_FT = max(0, NSIM_FT);
NSIM_MR = max(0, NSIM_MR);

%% Synthesize NSIM/STMI
% Regression of Fine-Timing NSIM + STMI performed by Rationalized Arcsine Transformation (RAU)
% Wirtzfeld et al (2017) 
% https://pmc.ncbi.nlm.nih.gov/articles/PMC5612921/

%SI_SYN = rau_b0 + rau_b1*STMI + rau_b2*NSIM_FT + rau_b3*STMI*NSIM_FT;

end