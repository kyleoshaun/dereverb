function [H] = MINT(G_rir, tau, results_dir)
% G_rir: Columns = acoustic impulse responses
% tau: Desired EIR delay (Impulse appears at index tau+1)

% L_h = AIR length (common to all channels, i.e., max AIR length, rest zero padded))
% M   = Number of channels (Microphones)
[L_h, M] = size(G_rir);

% Compute EQ FIR length (per channnel)
if (M > 1)
    L_g = ceil(((L_h-1) / (M-1)) * 1.1); % Invidual EQ FIR length
else
    L_g = L_h - 1;
end


% Construct multichannel sylvester filter matrix H = [H1 H2 ... HM]
G = [];
for ch_idx = 1:M
    g_rir_i = G_rir(:, ch_idx);
    G_i = sylvester_matrix(g_rir_i, L_g);

    G = [G G_i];
end

% Desired Equalized impulse response (EIR)
L_d = L_h + L_g - 1; % EIR length
d = zeros(L_d, 1); % EIR
d(tau+1) = 1;

% Compute MINT
h_mint = pinv(G) * d;

%fprintf("cond(H) = %d\n", cond(H))

%% Evaluate Equalized Impulse Response (EIR)

% Extract individual EQ Filters (per channel) and compute EIR
eir = zeros(L_d,1);
H = [];
for ch_idx = 1:M
    g_rir_i = G_rir(:, ch_idx);
    h_mint_i = h_mint((((ch_idx-1)*L_g)+1):(ch_idx*L_g));
    y_i = conv(g_rir_i, h_mint_i);

    H = [H h_mint_i];

    eir = eir + y_i;
end

g_rir_1 = G_rir(:, 1);
g_rir_2 = G_rir(:, 2);

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
subplot(2,1,1)
plot(g_rir_1)
title('RIR 1')
subplot(2,1,2)
plot(eir)
title('EIR')
sgtitle('MINT Results')

saveas(gcf, sprintf('%s/MINT.fig', results_dir));


