function MINT(h_air_list, tau)
% h_air_list: Columns = acoustic impulse responses
% tau: Desired EIR delay (Impulse appears at index tau+1)

% Shuffle and save seed for reproducibility
rng('shuffle')
seed = abs(round(randn(1,1) * 10000));
rng(seed);
%rng(724); % Reproduce specific seed

% L_h = AIR length (common to all channels, i.e., max AIR length, rest zero padded))
% M   = Number of channels (Microphones)
[L_h, M] = size(h_air_list);

% Compute EQ FIR length (per channnel)
if (M > 1)
    L_g = ceil(((L_h-1) / (M-1)) * 1.1); % Invidual EQ FIR length
else
    L_g = L_h - 1;
end


% Construct multichannel sylvester filter matrix H = [H1 H2 ... HM]
H = [];
for ch_idx = 1:M
    h_air_i = h_air_list(:, ch_idx);
    H_i = sylvester_matrix(h_air_i, L_g);

    H = [H H_i];
end

% Desired Equalized impulse response (EIR)
L_d = L_h + L_g - 1; % EIR length
d = zeros(L_d, 1); % EIR
d(tau+1) = 1;

% Compute MINT
g_mint = pinv(H) * d;

%fprintf("cond(H) = %d\n", cond(H))

%% Evaluate Equalized Impulse Response (EIR)

% Extract individual EQ Filters (per channel) and compute EIR
eir = zeros(L_d,1);
for ch_idx = 1:M
    h_air_i = h_air_list(:, ch_idx);
    g_mint_i = g_mint((((ch_idx-1)*L_g)+1):(ch_idx*L_g));
    y_i = conv(h_air_i, g_mint_i);

    eir = eir + y_i;
end

h_air_1 = h_air_list(:, 1);
h_air_2 = h_air_list(:, 2);

% figure()
% subplot(3,1,1)
% stem(h_air_1)
% title('AIR 1')
% subplot(3,1,2)
% stem(h_air_2)
% title('AIR 2')
% subplot(3,1,3)
% stem(eir)
% title('EIR')
% sgtitle('MINT Results')


figure()
subplot(3,1,1)
plot(h_air_1)
title('AIR 1')
subplot(3,1,2)
plot(h_air_2)
title('AIR 2')
subplot(3,1,3)
plot(eir)
title('EIR')
sgtitle('MINT Results')

