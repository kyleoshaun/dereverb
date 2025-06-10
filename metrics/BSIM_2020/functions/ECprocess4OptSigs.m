% Usage:[optsigmin,optsigmax] = ECprocess4OptSigs(optsigL,optsigR,ECparams,fs,bin_inaccuracy)
% This function applies the Equalization-Cancellation process in the frequency
% domain to the optional signals.
% The Equalization is applied two both ears symmetrically.
%
% Input:
% optsigL - left ear optional signal
% optsigR - right ear optional signal
% ECparams- Parameters for the EC process
% fs - sampling frequency
% bin_inaccuracy - flag indicating the use of (1) binaural processing
% inaccuracies or (0) assuming binaural processing to be deterministic
%
% Output: 
% optsigmin - EC processed optional signal using minimization strategy
% optsigmax - EC processed optional signal using maximization strategy
% Authors: Christopher Hauth <christopher.hauth@uni-oldenburg.de>
%          Dr. Thomas Brand  <thomas.brand@uni-oldenburg.de>
% Date: 24.05 2020
% Version: 1.0.
%--------------------------------------------------------------------------
function [optsigmin,optsigmax] = ECprocess4OptSigs(optsigL,optsigR,ECparams,fs,bin_inaccuracy)
orig_len = length(optsigL);
optsigLFFT = fft(optsigL,ECparams.fftpoints); 
optsigRFFT = fft(optsigR,ECparams.fftpoints);

f = linspace(0,fs/2,floor(ECparams.fftpoints./2+1))';
 % invert frequency vector for complex conjugate
f_inv = f(end-1:-1:2); 
 % frequency vector in rad
omega = 2.*pi*[f;f_inv];
 
 % copy signals to generate buffer for processed signals
FFToptsigLeq = optsigLFFT;
FFToptsigReq = optsigRFFT;
%--------------------------------------------------------------------------
                %% Apply Equalization Mechanism %%
% calculate phase factor for left and right ear with a set of gaussian
% distributed RVs as binaural processing inaccuracies
if bin_inaccuracy
    for kk = 1:size(ECparams.errorLR,2)
        phaseL(:,kk) = exp(-1j.*omega(1:(end/2)+1).*(ECparams.EC_tau./2 + ECparams.errorLR(1,kk)));
        phaseR(:,kk) = exp(+1j.*omega(1:(end/2)+1).*(ECparams.EC_tau./2 + ECparams.errorLR(2,kk)));
        % complex conjugate
        phaseLcc(:,kk) = exp(+1j.*omega((end/2)+2:end).*(ECparams.EC_tau./2 + ECparams.errorLR(1,kk)));
        phaseRcc(:,kk) = exp(-1j.*omega((end/2)+2:end).*(ECparams.EC_tau./2 + ECparams.errorLR(2,kk)));
    end
    % take the mean over all phase terms to obtain averaged processing
    % binaural processing inaccuracy
    phaseL   = mean(phaseL,2);
    phaseLcc = mean(phaseLcc,2);
    phaseR   = mean(phaseR,2);
    phaseRcc = mean(phaseRcc,2);
else
        phaseL = exp(-1j.*omega(1:(end/2)+1).*(ECparams.EC_tau./2));
        phaseR = exp(+1j.*omega(1:(end/2)+1).*(ECparams.EC_tau./2));
        % complex conjugate
        phaseLcc = exp(+1j.*omega((end/2)+2:end).*(ECparams.EC_tau./2));
        phaseRcc = exp(-1j.*omega((end/2)+2:end).*(ECparams.EC_tau./2));
end
% modify signals
% Left ear
FFToptsigLeq(1:(end/2)+1)  = FFToptsigLeq(1:(end/2)+1).*phaseL;
FFToptsigLeq((end/2)+2:end)= FFToptsigLeq((end/2)+2:end).*phaseLcc;
 
% Right ear 
FFToptsigReq(1:(end/2)+1)  = FFToptsigReq(1:(end/2)+1).*phaseR;
FFToptsigReq((end/2)+2:end)= FFToptsigReq((end/2)+2:end).*phaseRcc;
%--------------------------------------------------------------------------
                    %% Apply Cancellation %%
FFToptsigMin = FFToptsigLeq - FFToptsigReq; % destructive interference
FFToptsigMax = FFToptsigLeq + FFToptsigReq; % constructive interference
%--------------------------------------------------------------------------
                    %% Apply IFFT and cut signals
 % Level minimization
 optsigmin = real(ifft(FFToptsigMin));
 optsigmin = optsigmin(1:orig_len);
 
 % Level maximization
 optsigmax = real(ifft(FFToptsigMax));
 optsigmax = optsigmax(1:orig_len);              
%--------------------Licence ---------------------------------------------
% Copyright (c) <2020> Christopher Hauth
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