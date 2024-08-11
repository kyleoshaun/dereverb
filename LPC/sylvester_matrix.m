function H = sylvester_matrix(h, Lx)
% h = system IR
% L = Input signal length (to be filtered by matrix)

if isrow(h)
    h = h';
end

Lh = length(h);
Ly = Lh + Lx - 1;

H = zeros(Ly, Lx);
for n = 1:Lx
    H(n:(n+Lh-1),n) = h;
end
    
end