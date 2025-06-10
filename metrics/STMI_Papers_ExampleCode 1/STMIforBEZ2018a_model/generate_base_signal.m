function [ base_signal ] = generate_base_signal( sig )
%GENERATE_BASE_SIGNAL Generate a "base-signal" for the STMI calculation by
% randomizing the phase of a signal while keeping its long-term spectrum

D = size(sig);

if prod(D)>length(sig)
    error('Input signal must be a vector (i.e., one-dimensional array)')
end

Xk = fft(sig);

Xk = Xk(:);

n = length(sig);

if mod(n,2)==0
    f = [1:(n/2-1)]';
    p = 2*pi*rand(size(f));
    yc = Xk(2:(n/2)).*exp(1i.*p.*f);
    s = [Xk(1); yc; Xk(n/2+1); conj(flipud(yc))];
else
    f = [1:(n-1)/2]';
    p = 2*pi*rand(size(f));
    yc = Xk(2:((n-1)/2+1)).*exp(1i.*p.*f);
    s = [Xk(1); yc; conj(flipud(yc))];
end

x = ifft(s);
base_signal = real(x);

if D(2)>1
    base_signal = base_signal.';
end

end

