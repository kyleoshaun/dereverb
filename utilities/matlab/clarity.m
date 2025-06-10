function  C = clarity(rir, fs)
% Assumption: measurement noise has been removed (first non-zero RIR value is direct sound)

T_ER = 50 / 1000; % Duration of early reflections
N_ER = T_ER * fs;

% Locate direct sound
N_d = min(find(abs(rir) > 0));

if length(rir) > N_ER
    rir_useful = rir(N_d:(N_d+N_ER-1));
    rir_detrim = rir((N_d+N_ER):end);
    
    C = (sum(rir_useful .* rir_useful)) / (sum(rir_detrim .* rir_detrim));
else
    C = 10 ^ (40 / 20); % Saturate at 40 dB (C50 goes to infinity)
end

end