function [s_est, dap_equalizer_struct] = delay_and_predict(Y, s_ref, p1, p2, fs, s1_enable, s1_on_clean_speech, results_dir)
% Run Delay-and-Predict dereverberation with parallelized for loops
%
% Syntax
%   [s_est, dap_equalizer_struct] = delay_and_predict(Y, s_ref, p1, p2, fs, s1_enable, s1_on_clean_speech, results_dir)
%
% Inputs:
% - Y: Reverberant microphone signal matrix (Each column is a different microphone signal)
% - s_ref: Clean speech (only used if stage_1_on_clean_speech is enabled)
% - p1: Prediction order for source whitening stage
% - p2: Prediction order for multichannel linear prediction stage
% - fs: Sample rate (Hz)
% - s1_enable: Set to true to enable stage 1 of delay-and-predict (source whitening)
% - s1_on_clean_speech: Set to true to compute source whitening filter using clean speech (not blind)
% - results_dir: Directory to save results
%
% Outputs:
% - s_est: Estimate of clean speech (Blind etimate of s_ref from Y)
% - dap_equalizer_struct.H: Equalizer filters (Each column is the filter applied to the corresponding microphone signal)
% - dap_equalizer_struct.delays: Sample delays used per channel for time-alignment
% - dap_equalizer_struct.h0: Weights used for linear combination of individual MC-LP results

%% Create results folder

% Results can be saved in the Results Folder
if ~exist(results_dir,'Dir')
    mkdir(results_dir)
end


%% Parameters

M = size(Y,2); % Number of channels = number of columns (number of reverberant signals)

% Script Modes
apply_time_alignment    = false;

% Regularization Factor
% Diagonal reg_matrix=reg_factor*I is applied to the autocorrelation matrices
% used in both linear-prediction stages (i.e., R_ym and R_mc)
% to improve numerical stability of algorithm at the cost of reduced cancellation
% of lower energy reverb. This trade-off may be worth while since the low energy
% part of the reverb tail is often close in magnitude to the noise floor
% and therefore is cancelled effectively by the algorithm (in fact the resulting
% equalizer often introduces low energy reverb)
reg_factor = 5*10^-5; % T60_max = 1 sec
reg_factor = 2.5*10^-5; % T60_max = 500 msec
%reg_factor = 0;
reg_factor_1 = reg_factor;
reg_factor_2 = reg_factor;

% FFT/PSD params for plots
Nfft = 4096;
Nfft_welch  = 2^(nextpow2(length(s_ref)));
Nfft_welch  = min(Nfft_welch, Nfft);
Nwin_welch  = Nfft_welch / 2;
Nover_welch = Nfft_welch / 4;


%% Print Log
fprintf("\nDelay-and-Predict config:\n")
fprintf(" - Number of Microphones (M) = %.0f\n", M)
fprintf(" - Source whitening order (p1) = %.0f\n", p1)
fprintf(" - Multichannel Linear Prediction order (p2) = %.0f\n", p2)
fprintf(" - Source whitening Enabled? = %.0f\n", s1_enable)
fprintf(" - Source whitening on clean speech? = %.0f\n", s1_on_clean_speech)
fprintf(" - Regularization Factor for Source Whitening = %d\n", reg_factor_1)
fprintf(" - Regularization Factor for MC-LP = %d\n", reg_factor_2)
fprintf(" - Time Alignment Applied? = %.0f\n\n", apply_time_alignment)

%% LPC on clean speech (Not actually used in algo, just for reference)

alpha_s = 0; % Placeholder for when this is disabled

if s1_on_clean_speech

fprintf("\nRunning Stage 1 (Source Whitening) on CLEAN SPEECH ... ");


L          = length(s_ref);
w_hamming  = hamming(L);
s_w = s_ref .* w_hamming;

% Zero padding (for autocorr method) -- not sure this is required but doesnt hurt
s_w = [zeros(p1,1) ; s_w ; zeros(p1,1)];

% Biased autocorr calc because unbiased xcorr introduces a changing scale factor (n-|m|)
% which differs from the set up of the autocorrelation method for LPC
% (required for a stable inverse filter)
[phi_s, lags] = xcorr(s_w, s_w, 'biased', p1); 
idx_lag0    = find(lags==0);

R_s = toeplitz(phi_s(idx_lag0:(idx_lag0+p1-1)));
r_s = phi_s((idx_lag0+1):(idx_lag0+p1));

% Apply regularization factor to autocorrelation matrix
reg_matrix = reg_factor_1 .* eye(size(R_s));
R_s       = R_s + reg_matrix;

alpha_s = R_s \ r_s; % Option 1: Solve by gaussian elimination
%alpha_s = pinv(R_s) * r_s; % Option 2: Solve by pseudo-inverse


% Compute source-whitened reverberant microphone signals
X = filter([1 ; -1*alpha_s], 1, Y); % Apply source-whitening prediction error filter to each column

% Apply prediction error filter to clean speech to see how well source is whitened (for ref only)
s_ref_whitened = filter([1 ; -1*alpha_s], 1, s_ref);

fprintf("Done")

% Clear large variables no longer needed
varData = whos('R_s');
memory_info.R_s_size  = size(R_s);
memory_info.R_s_bytes = varData.bytes;
memory_info.R_s_type  = varData.class;
memory_info.R_s_cond  = cond(R_s);
clear R_s
clear r_s
clear phi_s

%% Stage 1: Source Whitening

alpha_ym = 0; % Placeholder for when this is disabled

else

fprintf("\nRunning Stage 1 (Source Whitening) ... ");

L          = length(s_ref);
w_hamming  = hamming(L);

% Apply hamming window to each signal, and zero pad (autocorrelation method)
Y_w = zeros(size(Y) + [2*p1 0]);
for ch = 1:M
    Y_w(:,ch) = [zeros(p1,1) ; Y(:,ch) .* w_hamming; zeros(p1,1)];
end

% Compute autocorrelation for each reverberant microphone signal
phi_Y = zeros(2*p1+1, M);
for ch = 1:M
    [phi_Y(:,ch), lags] = xcorr(Y_w(:,ch), 'biased', p1);
end
idx_lag0     = find(lags==0);

phi_ym = sum(phi_Y, 2); % Avg: Sum all autocorrelation functions (sum columns)

R_ym = toeplitz(phi_ym(idx_lag0:(idx_lag0+p1-1)));
r_ym = phi_ym((idx_lag0+1):(idx_lag0+p1));

% Apply regularization factor to autocorrelation matrix
reg_matrix = reg_factor_1 .* eye(size(R_ym));
R_ym       = R_ym + reg_matrix;

alpha_ym = R_ym \ r_ym; % Option 1: Solve by Gaussian elimination
%alpha_ym = pinv(R_ym) * r_ym; % Option 2: Solve by pseudo-inverse

%e_y1_m = filter([1 ; -1*alpha_ym], 1, y1_frame);

% Compute source-whitened reverberant microphone signals
X = filter([1 ; -1*alpha_ym], 1, Y);  % Apply source-whitening prediction error filter to each column

% Apply prediction error filter to clean speech to see how well source is whitened (for ref only)
s_ref_whitened = filter([1 ; -1*alpha_ym], 1, s_ref);

fprintf("Done")

% Clear large variables no longer needed
varData = whos('R_ym');
memory_info.R_ym_size  = size(R_ym);
memory_info.R_ym_bytes = varData.bytes;
memory_info.R_ym_type  = varData.class;
memory_info.R_ym_cond  = cond(R_ym);
clear R_ym
clear r_ym
clear phi_ym

end

%% Stage 1 Plots (For first 2 channels)

y1 = Y(:,1);
y2 = Y(:,2);
x1 = X(:,1);
x2 = X(:,2);

[Sm, freqs_welch] = pwelch(s_ref, Nwin_welch, Nover_welch, Nfft_welch, fs);
Y1m = pwelch(y1, Nwin_welch, Nover_welch, Nfft_welch);
Y1m = pwelch(y2, Nwin_welch, Nover_welch, Nfft_welch);
X1m = pwelch(x1, Nwin_welch, Nover_welch, Nfft_welch);
X2m = pwelch(x2, Nwin_welch, Nover_welch, Nfft_welch);
S_ref_whitened_m = pwelch(s_ref_whitened, Nwin_welch, Nover_welch, Nfft_welch);


if s1_on_clean_speech
    [V_s, w_freqz] = freqz(1, [1 ; -1*alpha_s], Nfft); 
    freqs_freqz   = w_freqz * (fs / (2*pi));

    figure()
    subplot(2,1,1)
    plot(freqs_welch ./ 1000, 10*log10(Sm));
    hold on;
    plot(freqs_freqz ./ 1000,20*log10(abs(V_s .* (sqrt(max(Sm)) / max(abs(V_s))))));
    xlabel('Frequency [kHz]')
    ylabel('dB')
    %title('Vocal Tract model')
    legend('Clean Speech Spectrum (Not reverberant)', 'LPC Inverse Filter')
    title('Source Spectrum Estimation')
    subplot(2,1,2)
    plot(freqs_welch ./ 1000, 10*log10(S_ref_whitened_m));
    xlabel('Frequency [kHz]')
    ylabel('dB')
    title('Whitened Source Spectrum')
    sgtitle("Stage 1: Results of LPC on Clean Speech (NOT BLIND)")

    saveas(gcf, sprintf('%s/S1_SourceWhiteningLPC.fig', results_dir));
else
    [V_ym, w_freqz] = freqz(1, [1 ; -1*alpha_ym], Nfft); 
    freqs_freqz   = w_freqz * (fs / (2*pi));

    figure()
    subplot(2,1,1)
    plot(freqs_welch ./ 1000, 10*log10(Sm));
    hold on;
    plot(freqs_freqz ./ 1000,20*log10(abs(V_ym .* (sqrt(max(Sm)) / max(abs(V_ym))))));
    xlabel('Frequency [kHz]')
    ylabel('dB')
    %title('Vocal Tract model')
    legend('Clean Speech Spectrum (Not reverberant)', 'LPC Inverse Filter')
    title('Spectrum Estimation')
    subplot(2,1,2)
    plot(freqs_welch ./ 1000, 10*log10(S_ref_whitened_m));
    xlabel('Frequency [kHz]')
    ylabel('dB')
    title('Whitened Clean Speech Spectrum')
    sgtitle("Stage 1: Results of LPC on Reverberant Speech (Blind, Avg over channels)")

    saveas(gcf, sprintf('%s/S1_SourceWhiteningLPC.fig', results_dir));

end

% %% Analysis of conditioning at every stage
% 
% y1_w        = y1_frame .* w_hamming;
% x1_w        = x1_frame .* w_hamming;
% 
% y1_w = [zeros(p1,1) ; y1_w ; zeros(p1,1)];
% x1_w = [zeros(p1,1) ; x1_w ; zeros(p1,1)];
% 
% [phi_y1y1, lags] = xcorr(y1_w, y1_w, 'biased', p1);
% [phi_x1x1,    ~] = xcorr(x1_w, x1_w, 'biased', p1);
% idx_lag0         = find(lags==0);
% 
% R_y1 = toeplitz(phi_y1y1(idx_lag0:(idx_lag0+p1-1)));
% R_x1 = toeplitz(phi_x1x1(idx_lag0:(idx_lag0+p1-1)));
% 
% figure()
% plot(freqs_welch ./ 1000, 10*log10(Sm_frame));
% xlabel('Frequency [kHz]')
% ylabel('dB')
% title('Clean Spectrum')
% 
% figure()
% plot(freqs_welch ./ 1000, 10*log10(Esm_frame));
% xlabel('Frequency [kHz]')
% ylabel('dB')
% title('Whitened Clean Spectrum')
% 
% figure()
% plot(freqs_welch ./ 1000, 10*log10(Y1m_frame));
% xlabel('Frequency [kHz]')
% ylabel('dB')
% title('Reverberant Spectrum')
% 
% figure()
% plot(freqs_welch ./ 1000, 10*log10(X1));
% xlabel('Frequency [kHz]')
% ylabel('dB')
% title('Source-Whitened Reverberant Spectrum')
% 
% fprintf("Condition Numbers:\n");
% fprintf("cond(R_s) = %d\n", cond(R_s));
% fprintf("cond(R_y) = %d\n", cond(R_y1));
% fprintf("cond(R_x) = %d\n", cond(R_x1));

%% STAGE 2: MULTI-CHANNEL LINEAR PREDICTION ON SOURCE-WHITENED REVERBERANT SIGNALS

fprintf("\nRunning Stage 2 (Multichannel Linear Prediction) ... ")

if s1_enable == false
    % Perform MC-LPC directly on reverberant signals (not source-whitened)
    X = Y;
end

h0 = zeros(M,1); % First vector coefficient (estimate) of multi-path channel -- for time alignment
delay_comp = -1;
delays = zeros(M,1);
%delay_filters = ;
%delay_filter_1 = [1];
%delay_filter_2 = [1];
%delay_filter_3 = [1];
%delay_filter_4 = [1];
%delay_filter_5 = [1];
%delay_filter_6 = [1];
coeff_thresh = 0.1;

while min(abs(h0)) < coeff_thresh

    delay_comp = delay_comp + 1;

    %fprintf("\n  Iteration %.0f ... ", delay_comp);

    for ch = 1:M
        if abs(h0(ch)) > coeff_thresh
            X(:,ch) = [0 ; X(1:(end-1), ch)]; % Delay signal by one sample
            delays(ch) = delays(ch) + 1;
        end
    end

    L          = length(X(:,1));
    w_hamming  = hamming(L);


    % Apply hamming window to each signal, and zero pad (autocorrelation method)
    X_w = zeros(size(X) + [2*p2 0]);
    for ch = 1:M
        X_w(:,ch) = [zeros(p2,1) ; X(:,ch) .* w_hamming; zeros(p2,1)];
    end

    % Compute all auto/cross correlations (every combination of signals)
    [phi_X, lags] = xcorr(X_w, 'biased', p2);
    idx_lag0       = find(lags==0);
    % I.e.,  phi_X = [phi_x1x1  phi_x1x2 ... phi_x1xM ... phi_xMxM]
    %                 :  |         |            |            |
    %               lags |         |            |            |
    %                 :  |         |            |            |
    %
    % Note: Rows = lags, col = M * (first channel idx - 1) + (second channel idx)

    % Initialize variables for full multi-channel spatio-temporal correlation matrix / block-vector
    R_mc = zeros(M*p2, M*p2);
    r_mc = zeros(M,    M*p2);
    % R_mc = [ R_xx(0)       R_xx(1)     ...  R_xx(p-1) ]
    %        [ R_xx(-1)      R_xx(0)     ...  R_xx(p-2) ]
    %        [    :             :         \      :      ]
    %        [ R_xx(-p2+1)  R_xx(-p2+1)  ...  R_xx(0)   ]
    % r_mc = [ R_xx(1)  R_xx(2)  ...  R_xx(p) ]
    % Note: Although r_mc is a matrix, it has a single row of size-MxM spatial correlation matrices (blocks)
    
    
    % Construct full multi-channel spatio-temporal correlation matrix (R_mc)
    %k = 0; 
    ks = 0:-1:(-p2+1);
    R_mc_block_rows = cell(p2,1);
    parfor block_row = 1:p2
        k = ks(block_row);
        R_mc_block_row = zeros(M, M*p2);
        for block_col = 1:p2
            R_xx_k = zeros(M,M); % Spatial (inter-mic) cross-correlation matrix for a lag of k
            for ch1 = 1:M
                for ch2 = 1:M
                    phi_X_col = M*(ch1-1) + ch2;
                    R_xx_k(ch1, ch2) = phi_X(idx_lag0+k, phi_X_col);
                end
            end

            % E.g., 6 channel: R_xx(k) (denoted R_xx_k)
            % R_xx_k = [phi_x1x1(idx_lag0+k) phi_x1x2(idx_lag0+k) phi_x1x3(idx_lag0+k) phi_x1x4(idx_lag0+k) phi_x1x5(idx_lag0+k) phi_x1x6(idx_lag0+k) ; ...
            %           phi_x2x1(idx_lag0+k) phi_x2x2(idx_lag0+k) phi_x2x3(idx_lag0+k) phi_x2x4(idx_lag0+k) phi_x2x5(idx_lag0+k) phi_x2x6(idx_lag0+k) ; ...
            %           phi_x3x1(idx_lag0+k) phi_x3x2(idx_lag0+k) phi_x3x3(idx_lag0+k) phi_x3x4(idx_lag0+k) phi_x3x5(idx_lag0+k) phi_x3x6(idx_lag0+k) ; ...
            %           phi_x4x1(idx_lag0+k) phi_x4x2(idx_lag0+k) phi_x4x3(idx_lag0+k) phi_x4x4(idx_lag0+k) phi_x4x5(idx_lag0+k) phi_x4x6(idx_lag0+k) ; ...
            %           phi_x5x1(idx_lag0+k) phi_x5x2(idx_lag0+k) phi_x5x3(idx_lag0+k) phi_x5x4(idx_lag0+k) phi_x5x5(idx_lag0+k) phi_x5x6(idx_lag0+k) ; ...
            %           phi_x6x1(idx_lag0+k) phi_x6x2(idx_lag0+k) phi_x6x3(idx_lag0+k) phi_x6x4(idx_lag0+k) phi_x6x5(idx_lag0+k) phi_x6x6(idx_lag0+k)];
    
            col_0 = (block_col - 1) * M + 1;
            col_1 = (block_col) * M;

            R_mc_block_row(:, (col_0:col_1)) = R_xx_k;

            k = k + 1;
        end
        R_mc_block_rows{block_row} = R_mc_block_row;
    end

    % Combine block rows into full multi-channel spatio-temporal correlation matrix (R_mc)
    for block_row = 1:p2
        row_0 = (block_row - 1) * M + 1;
        row_1 = (block_row) * M;
        
        R_mc((row_0:row_1), :) = R_mc_block_rows{block_row};
    end

    % Construct full multi-channel spatio-temporal correlation block-vector (r_mc)
    k = 1;
    for block_col = 1:p2

        R_xx_k = zeros(M,M); % Spatial (inter-mic) cross-correlation matrix for a lag of k (Same as before)
        for ch1 = 1:M
            for ch2 = 1:M
                phi_X_col = M*(ch1-1) + ch2;
                R_xx_k(ch1, ch2) = phi_X(idx_lag0+k, phi_X_col);
            end
        end  

        % Insert into full multi-channel spatio-temporal correlation block-vector (r_mc)
        col_0 = (block_col - 1) * M + 1;
        col_1 = (block_col) * M;
    
        r_mc(:, col_0:col_1) = R_xx_k;
        k = k+1;
    end

    % Apply regularization factor to autocorrelation matrix
    reg_matrix = reg_factor_2 .* eye(size(R_mc));
    R_mc       = R_mc + reg_matrix;
    
    % Solve Multichannel Normal Equations
    % Note row formulation of system is used here, so aR = r (rather than Ra = r)
    %alpha_mc = r_mc / R_mc; % Option: Solve by Gaussian Elimination
    alpha_mc = r_mc * pinv(R_mc); % Option 2: Solve by pseudo inverse
    %alpha_mc = r_mc * ((inv(R_mc'*R_mc))*R_mc'); % Pseudo inverse for over-determined system
    %alpha_mc = r_mc * (R_mc'*(inv(R_mc*R_mc'))); % Pseudo inverse for under-determined system
  

    % Clear large variables no longer needed
    varData = whos('R_mc');
    memory_info.R_mc_size  = size(R_mc);
    memory_info.R_mc_bytes = varData.bytes;
    memory_info.R_mc_type  = varData.class;
    memory_info.R_mc_cond   = cond(R_mc);
    clear R_mc
    clear r_mc
    clear phi_X

    % Extract individual (single-channel) prediction Filters
    % and put them into the rows of a matrix H
    H = zeros(p2+1, M*M);

    % h_pred_i_from_k = [0 alpha_mc(i, k:M:((p2-1)*M+k))];
    for ch_i = 1:M
        for ch_k = 1:M
            H_col = M*(ch_i-1) + ch_k;
            h_pred_i_from_k = [0 alpha_mc(ch_i, ch_k:M:((p2-1)*M+ch_k))]';
            H(:,H_col) = h_pred_i_from_k;
        end
    end
    % I.e.,  H = [h_pred_1_from_1  h_pred_1_from_2 ... h_pred_1_from_M  h_pred_2_from_1  ...  h_pred_M_from_M]
    %              :    |                 |                  |                 |                     |
    %      filter taps  |                 |                  |                 |                     |
    %              :    |                 |                  |                 |                     |
    %
    % Note: Rows = filter taps, col = M * (channel being predicted - 1) + (Channel used to predict)    

    % Compute/combine individual (single-channel) prediction filter outputs to form multichannel prediction filter
    % I.e., xi_est = h_pred_i_from_1 * x1 + h_pred_i_from_2 * x2 +  ... + h_pred_i_from_M * xM
    X_est = zeros(size(X));
    for ch_i = 1:M
        xi_est = zeros(length(X(:,1)), 1);
        for ch_k = 1:M
            xk = X(:, ch_i);
            H_col = M*(ch_i-1) + ch_k;
            h_pred_i_from_k = H(:,H_col);
            xi_est = xi_est + filter(h_pred_i_from_k, 1, xk);
        end

        % Collect results of each prediction in a data matrix
        X_est(:, ch_i) = xi_est;
    end

    % Compute Corresponding Prediction errors
    % I.e., e_xi_est = xi - xi_est
    E_est = X - X_est;
    
    % Estimate first vector coefficient of SIMO channel (for linear combiner)
    L          = length(E_est(:,1));
    w_hamming  = hamming(L);
    X_res_w = zeros(size(X)); % res: Residual
    for ch = 1:M
        X_res_w(:,ch) = E_est(:,ch) .* w_hamming; % Apply hamming window to each column
    end

    % Compute all auto/cross correlations (every combination of signals)
    [phi_res, lags] = xcorr(X_res_w, 'biased', p2);
    idx_lag0       = find(lags==0);
    % Note: Same structure as before (see phi_X)
    
    % Extract spatial correlation matrix of multichannel residual (prediction error) for lag k=0
    k = 0;
    R_res_0 = zeros(M,M);
    for ch1 = 1:M
        for ch2 = 1:M
            phi_res_col = M*(ch1-1) + ch2;
            R_res_0(ch1, ch2) = phi_res(idx_lag0+k, phi_res_col);
        end
    end
    
    [V,D] = eig(R_res_0);
    h0 = V(:, find(diag(D) == max(diag(D)))); % Estimated vector coefficient of SIMO channel

    %fprintf("Done");

    if apply_time_alignment ~= true
        break;
    end
end

Y_delayed = Y;
if apply_time_alignment
    max_delay = max(delays);

    % delay_filter_k = [zeros(delays(k),1) ; 1 ; zeros((max_delay - delays(k)),1)];
    %for k = 1:M
    %    delay_filter_k = [zeros(delays(k),1) ; 1 ; zeros((max_delay - delays(k)),1)];
    %    delay_filters(:,k) = delay_filter_k;
    %end

    % Adjust other signals for delay compensation
    s_ref  = filter([zeros(delay_comp,1) ; 1], 1, s_ref);
    for k = 1:M
        %delay_filter_k = delay_filters(:,k)
        %Y_delayed(:,k) = filter(delay_filter_k, 1, Y_delayed(:,k));
        Y_delayed(:,k) = filter([zeros(delays(k),1) ; 1], 1, Y_delayed(:,k))
    end
end

fprintf("Done")

%% Create equalizer structure (for output)
dap_equalizer_struct.H          = H;
dap_equalizer_struct.delays     = delays;
dap_equalizer_struct.h0         = h0;
dap_equalizer_struct.delay_comp = delay_comp;


%% Stage 3: Inverse Filtering

fprintf("\nRunning Stage 3 (Inverse Filtering) ... ")

% Compute DAP output (de-reverberate)

% Re-apply multi-channel prediction error filter to 
% original reverberant microphones signals (not source-whitened)
s_est = apply_dap_equalizer(dap_equalizer_struct, Y);

fprintf("Done\n\n")

%% Save results
save(sprintf('%s/X.mat', results_dir),'X');
save(sprintf('%s/X.mat', results_dir),'X');
save(sprintf('%s/Y.mat', results_dir),'Y');
save(sprintf('%s/dap_equalizer_struct.mat', results_dir),'dap_equalizer_struct');
%save(sprintf('%s/R_mc.mat', results_dir),'R_mc', '-v7.3'); % MAT v7.3 for large matrix
%save(sprintf('%s/r_mc.mat', results_dir),'r_mc');
%save(sprintf('%s/R_ym.mat', results_dir),'R_ym', '-v7.3'); % MAT v7.3 for large matrix
%save(sprintf('%s/r_ym.mat', results_dir),'r_ym');
%save(sprintf('%s/R_s.mat', results_dir),'R_s', '-v7.3'); % MAT v7.3 for large matrix
%save(sprintf('%s/r_s.mat', results_dir),'r_s');
save(sprintf('%s/memory_info.mat', results_dir),'memory_info');
save(sprintf('%s/alpha_s.mat', results_dir),'alpha_s');
save(sprintf('%s/alpha_ym.mat', results_dir),'alpha_ym');
    
end
