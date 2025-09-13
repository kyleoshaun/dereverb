function [s_est] = apply_dap_equalizer(dap_equalizer_struct, Y)
% Apply delay-and-predict equalizer to blindly dereverberate reverberant speech
% Inputs:
% - dap_equalizer_struct (Equalizer object returned by delay_and_predict()
% - Y: Reverberant microphone signal matrix (Each column is a different microphone signal)
% Outputs:
% - s_est: Dereverberated Speech (i.e., estimate of clean speech)

H      = dap_equalizer_struct.H;
delays = dap_equalizer_struct.delays;
h0     = dap_equalizer_struct.h0;

M = size(Y,2); % Number of channels = number of columns (number of reverberant signals)

% Apply Delays for time-alignment
Y_delayed = zeros(size(Y));
for k = 1:M
    Y_delayed(:,k) = filter([zeros(delays(k),1) ; 1], 1, Y(:, k));
end

% Compute/combine individual (single-channel) prediction filter outputs to form multichannel prediction filter
% I.e., yi_est = h_pred_i_from_1 * y1 + h_pred_i_from_2 * y2 +  ... + h_pred_i_from_M * yM
Y_est = zeros(size(Y));
for ch_i = 1:M
    yi_est = zeros(length(Y(:,1)), 1);
    for ch_k = 1:M
        yk = Y(:, ch_k);
        H_col = M*(ch_i-1) + ch_k;
        h_pred_i_from_k = H(:,H_col);
        yi_est = yi_est + filter(h_pred_i_from_k, 1, yk);
    end

    % Collect results of each prediction in a data matrix
    Y_est(:, ch_i) = yi_est;
end

% Compute Corresponding Prediction errors, which should represent the dereverberated
% signal (the clean source)
% I.e., si_est = si - si_est
S_est = Y - Y_est;

% Apply linear combiner between all M dereverberated signals
s_est = S_est * h0;

end