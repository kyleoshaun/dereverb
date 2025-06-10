function [edc] = EDC(rir)
total_energy_rir = rir' * rir;
remaining_energy_rir = total_energy_rir;
edc    = zeros(length(rir), 1);
for n = 1:length(rir)
    remaining_energy_rir = (rir(n:end))' * (rir(n:end));
    edc(n) = remaining_energy_rir / total_energy_rir;
end
end