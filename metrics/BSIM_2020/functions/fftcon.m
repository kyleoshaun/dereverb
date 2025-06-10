% USAGE: [Mixsigmin, Mixsigmax, ECparams4Opt] = fftcon(MixsigL,MixsigR,fc,fs,sigmadelta0,Delta0,bin_inaccuracy)
% Description:
% fftcon.m calculates the interaural delay of the dominant source in the
% frequency domain and applies the Cancellation in the subesequent function
% EC_FFTprocess.m. 
% First, left and right ear signal are transformed to
% frequency domain in order to calculate the cross power spectral density.
% By calculating the angle of the CPSD and, the frequency dependent ITD can
% be obtained by dividing the angle by the angular frequency. In order to
% obtain the final ITD value for the frequency band, an integration window
% derived from the absolute value of the CPSD is applied.
% Input: 
% MixsigL - Left ear signal
% MixsigR - Right ear signal
% fc      - center frequency of frequency channel
% fs      - sampling frequency
% sigmadelta0 - term to calculate binaural processing inaccuracy
% Delta0      - term to calculate binaural processing inaccuracy
% bin_inaccuracy - flag indicating the use of (1) binaural processing
% inaccuracies or (0) assuming binaural processing to be deterministic

% Output:
% Mixsigmin - EC processed signal using the minimization strategy
% Mixsigmax - EC processed signal using the maximization strategy
% ECparams4Opt - Structure containing EC parameters (Delay and Uncertainties)

% Author: Christopher F. Hauth <christopher.hauth@uni-oldenburg.de>
%         Dr. Thomas Brand     <thomas.brand@uni-oldenburg.de>
% Date: 24.05.2020
% Version: 1.0
function [Mixsigmin, Mixsigmax, ECparams4Opt] = fftcon(MixsigL,MixsigR,fc,fs,sigmadelta0,Delta0,bin_inaccuracy)

% caluclate number of fft points for computation of correlation
lenSig = length(MixsigL);
fftpoints = (length(MixsigL)+length(MixsigR));
% compute fft representation
MixsigLFFT = fft(MixsigL,fftpoints); 
MixsigRFFT = fft(MixsigR,fftpoints);

%% Frequency domain fractionaly delay estimation
  % Apply ERB window in the frequency domain
ERB       = 24.7.*(4.37.*fc./1000 + 1);
 % compute CrossPowerSpectralDensity
PHY_LR    = MixsigLFFT.*conj(MixsigRFFT);

Frequency = linspace(0, 44100, fftpoints);
PHY_phase = angle(PHY_LR);

% Make sure phase is consistent in case of pi (+-pi)
PHY_phase_unwrapped = PHY_phase;
PHY_phase_unwrapped(PHY_phase<0) = PHY_phase(PHY_phase<0)+2*pi;
% Convert IPD to ITD
PHY_time      = PHY_phase./(2.*pi.*Frequency');
PHY_time_unwrapped      = PHY_phase_unwrapped./(2.*pi.*Frequency');

% Normalize
PHY_abs_norm = abs(PHY_LR)./max(abs(PHY_LR));
window = find(Frequency>=fc-ERB/2&Frequency<=fc+ERB/2);

% compute standard deviation of tau 
EC_std = std(PHY_time(window));
% do the same for unwrapped phase:
EC_std_unwrapped = std(PHY_time_unwrapped(window));
% compute standard deviation of the delay processign inaccuracy
if  EC_std <= EC_std_unwrapped
    EC_tau = sum(PHY_time(window).*(PHY_abs_norm(window)/sum(PHY_abs_norm(window))));
else
    EC_tau = sum(PHY_time_unwrapped(window).*(PHY_abs_norm(window)/sum(PHY_abs_norm(window)))); %EC_tau_unwrapped 
end
std_min_frac = sigmadelta0.*(1+(abs(EC_tau)./(Delta0)));
% compute a set of 2 Gaussian variables for each ear
errorLR = 1.*std_min_frac.*randn(2,1);
% Save EC params in a struct. It will also be used to EC process optional
% signals
ECparams4Opt.errorLR = errorLR;
ECparams4Opt.fftpoints = fftpoints; 
ECparams4Opt.EC_tau = EC_tau; 
% apply EC mechanism in the frequency domain:
[Mixsigmin, Mixsigmax] =EC_FFTprocess(MixsigLFFT,MixsigRFFT,EC_tau,errorLR,fs,fftpoints,lenSig,bin_inaccuracy);
%--------------------Licence ---------------------------------------------
% Copyright (c) <2020> Christopher F. Hauth
% Dept. Medical Physics and Acoustics
% Carl von Ossietzky University Oldenburg 
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to 
% permit persons to whom the Software is furnished to do so, subject 
% to the following conditions:
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% END OF FILE