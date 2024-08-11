addpath('../../samples/')
addpath('../../AIR_Databases/AIR_1_4_BinauralDatabase/')

close all
clear
clc

% Shuffle and save seed for reproducibility
rng('shuffle')
seed = abs(round(randn(1,1) * 10000));
rng(seed);
%rng(724); % Reproduce specific seed

%% Parameters

M   = 4; % Num mics -- Note need to change the rest of the code to match this if changing
L_h = 1000; % Invidual Channel FIR length
if (M > 1)
    L_g = (L_h-1) / (M-1); % Invidual EQ FIR length
else
    L_g = L_h - 1;
end


%% Channel

% Exponential decay curve
tau = L_h/4;
exp_decay = exp(-1 .* (1:L_h)' ./ tau);

% Option 1: Generate random AIR
b_air_1 = randn(L_h,1);
b_air_2 = randn(L_h,1);
b_air_3 = randn(L_h,1);
b_air_4 = randn(L_h,1);

% Apply synthetic exponential decay to AIRs
b_air_1 = b_air_1(1:L_h) .* exp_decay;
b_air_2 = b_air_2(1:L_h) .* exp_decay;
b_air_3 = b_air_3(1:L_h) .* exp_decay;
b_air_4 = b_air_4(1:L_h) .* exp_decay;

% Synthetic 2-sided AIR
%b_air_1 = [flip(b_air_1) ; b_air_1];
%b_air_2 = [flip(b_air_2) ; b_air_2];

% % Add Synthetic time delay
% delay_1 = 0;
% delay_2 = 0;
% b_air_1 = [zeros(delay_1,1) ; b_air_1 ; zeros(delay_2,1)];
% b_air_2 = [zeros(delay_2,1) ; b_air_2  ; zeros(delay_1,1)];

% Option 2: Real AIR
% [b_air_1, b_air_2] = load_real_AIR(fs);
% a_air_1 = 1;
% a_air_2 = 1;

%% MINT

h_air_1 = b_air_1;
h_air_2 = b_air_2;
h_air_3 = b_air_3;
h_air_4 = b_air_4;

L_d = L_h + L_g - 1;

% Construct Sylvester Filter Matrix
H1 = sylvester_matrix(h_air_1, L_g);
H2 = sylvester_matrix(h_air_2, L_g);
H3 = sylvester_matrix(h_air_3, L_g);
H4 = sylvester_matrix(h_air_4, L_g);

H  = [H1 H2 H3 H4];

tau = 0;
d = zeros(L_d, 1);
d(tau+1) = 1;

g_mint = pinv(H) * d;

fprintf("cond(H) = %d\n", cond(H))

g1_mint = g_mint((((1-1)*L_g)+1):(1*L_g));
g2_mint = g_mint((((2-1)*L_g)+1):(2*L_g));
g3_mint = g_mint((((3-1)*L_g)+1):(3*L_g));
g4_mint = g_mint((((4-1)*L_g)+1):(4*L_g));

%% Evaluate Equalized Impulse Response (EIR)

y1 = filter(h_air_1, 1, g1_mint);
y2 = filter(h_air_2, 1, g2_mint);
y3 = filter(h_air_3, 1, g3_mint);
y4 = filter(h_air_4, 1, g4_mint);

eir = y1 + y2 + y3 + y4;

figure()
subplot(3,1,1)
stem(h_air_1)
title('AIR 1')
subplot(3,1,2)
stem(h_air_2)
title('AIR 2')
subplot(3,1,3)
stem(eir)
title('EIR')

