function [H] = MINT(G_rir, tau, fs, L_g, results_dir)
% G_rir: Columns = acoustic impulse responses
% tau: Desired EIR delay (Impulse appears at index tau+1)
% fs: sample rate (Hz)
% L_g: Per-channel FIR EQ length

% L_h = AIR length (common to all channels, i.e., max AIR length, rest zero padded))
% M   = Number of channels (Microphones)
[L_h, M] = size(G_rir);

% Compute EQ FIR length (per channnel)
% if (M > 1)
%     L_g = ceil(((L_h-1) / (M-1)) * 1.1); % Invidual EQ FIR length
% else
%     L_g = L_h - 1;
% end


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


ylim_max = 1.2*max(abs(g_rir_1));
figure()
subplot(2,1,1)
plot((0:(length(g_rir_1)-1)) .* (1/fs), g_rir_1)
ylim([-ylim_max ylim_max])
xlabel('Time [sec]')
title('Channel 1 Impulse Response')
subplot(2,1,2)
plot((0:(length(eir)-1)) .* (1/fs), eir)
ylim([-ylim_max ylim_max])
xlabel('Time [sec]')
title('Equalized Impulse Response (MINT)')

saveas(gcf, sprintf('%s/MINT.fig', results_dir));

edc_rir  = EDC(g_rir_1);
edc_MINT = EDC(eir);

figure()
plot((0:(length(edc_rir)-1)) .* (1/fs), 10*log10(edc_rir));
hold on;
plot((0:(length(edc_MINT)-1)) .* (1/fs), 10*log10(edc_MINT .* (max(edc_rir) / max(edc_MINT))));
ylim([-70 6])
%xlim([0 (length(h_channel_1) * (1/fs))])
xlabel('Time [sec]')
ylabel('Energy remaining [dB]')
legend('Reverb Energy', 'Equalized Reverb Energy')
title('Energy Decay Curve (MINT)')

saveas(gcf, sprintf('%s/EDC_MINT.fig', results_dir));

end