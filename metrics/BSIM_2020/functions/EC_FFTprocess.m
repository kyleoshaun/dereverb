% Usage: [Mixsigmin Mixsigmax] = EC_FFTprocess(MixsigLFFT,MixsigRFFT,EC_tau,errorLR,fs,fftpoints,orig_len,bin_inaccuracy)
% This function applies the Equalization-Cancellation process in the frequency
% domain.
% The Equalization is applied two both ears symmetrically.
% Input: 
% MixsigLFFT - FFT representation of left ear signal 
% MixsigRFFT - FFT representation of right ear signal
% EC_tau     - Delay for the EC process
% errorLR    - Pair of binaural processing inaccuracies
% fs         - sampling frequency
% fftpoints  - fftpoints in the fft
% orig_len   - original length of the signal
% bin_inaccuracy - flag indicating the use of (1) binaural processing
% inaccuracies or (0) assuming binaural processing to be deterministic

% Output: 
% Mixsigmin - EC processed signal applying minimization strategy
% Mixsigmax - EC processed signal applying maximization strategy
% Authors: Christopher Hauth <christopher.hauth@uni-oldenburg.de>
%          Dr. Thomas Brand  <thomas.brand@uni-oldenburg.de>
% Date: 22.10 2020
% Version: 1.0.
%--------------------------------------------------------------------------
function [Mixsigmin Mixsigmax] = EC_FFTprocess(MixsigLFFT,MixsigRFFT,EC_tau,errorLR,fs,fftpoints,orig_len,bin_inaccuracy)

 f = linspace(0,fs/2,floor(fftpoints./2+1))';
 % invert frequency vector for complex conjugate
 f_inv = f(end-1:-1:2); 
 % frequency vector in rad
 omega = 2.*pi*[f;f_inv];
 % copy signals to generate buffer for processed signals
 FFTMixsigLeq = MixsigLFFT;
 FFTMixsigReq = MixsigRFFT;
%--------------------------------------------------------------------------
                %% Apply Equalization Mechanism %%
% calculate phase factor for left and right ear with a set of gaussian
% distributed RVs as binaural processing inaccuracies
if bin_inaccuracy
    for kk = 1:size(errorLR,2)
        phaseL(:,kk) = exp(-1j.*omega(1:(end/2)+1).*(EC_tau./2 + errorLR(1,kk)));
        phaseR(:,kk) = exp(+1j.*omega(1:(end/2)+1).*(EC_tau./2 + errorLR(2,kk)));
        % complex conjugate
        phaseLcc(:,kk) = exp(+1j.*omega((end/2)+2:end).*(EC_tau./2 + errorLR(1,kk)));
        phaseRcc(:,kk) = exp(-1j.*omega((end/2)+2:end).*(EC_tau./2 + errorLR(2,kk)));
    end
    % take the mean over all phase terms to obtain averaged processing
    % binaural processing inaccuracy
    phaseL   = mean(phaseL,2);
    phaseLcc = mean(phaseLcc,2);
    phaseR   = mean(phaseR,2);
    phaseRcc = mean(phaseRcc,2);
else
        phaseL = exp(-1j.*omega(1:(end/2)+1).*(EC_tau./2));
        phaseR = exp(+1j.*omega(1:(end/2)+1).*(EC_tau./2));
        % complex conjugate
        phaseLcc = exp(+1j.*omega((end/2)+2:end).*(EC_tau./2));
        phaseRcc = exp(-1j.*omega((end/2)+2:end).*(EC_tau./2));
end
% equalize phases of left and right ear signals
% Left ear
FFTMixsigLeq(1:(end/2)+1)  = FFTMixsigLeq(1:(end/2)+1).*phaseL;
FFTMixsigLeq((end/2)+2:end)= FFTMixsigLeq((end/2)+2:end).*phaseLcc;
   
% Right ear 
FFTMixsigReq(1:(end/2)+1)  = FFTMixsigReq(1:(end/2)+1).*phaseR;
FFTMixsigReq((end/2)+2:end)= FFTMixsigReq((end/2)+2:end).*phaseRcc;
%--------------------------------------------------------------------------
                    %% Apply Cancellation %%
FFTMixsigMin = FFTMixsigLeq - FFTMixsigReq; % destructive interference
FFTMixsigMax = FFTMixsigLeq + FFTMixsigReq; % constructive interference
%--------------------------------------------------------------------------
                    %% Apply IFFT and cut signals
 % Level minimization
 Mixsigmin = real(ifft(FFTMixsigMin));
 Mixsigmin = Mixsigmin(1:orig_len); 
 % Level maximization
 Mixsigmax = real(ifft(FFTMixsigMax));
 Mixsigmax = Mixsigmax(1:orig_len);             
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