function c = symphase(c,a,varargin)
%SYMPHASE  Change Gabor coefficients to symmetric phase
%   Usage:  c=symphase(c,a);
%
%   SYMPHASE(c,a) alters the phase of the Gabor coefficients c so as if
%   they were obtained from a Gabor transform based on symmetric
%   time/frequency shifts. The coefficient must have been obtained from a
%   DGT with parameter a.
%
%   Gabor coefficients with symmetric phase correspond to the following 
%   transform:
%   Consider a signal f of length L and define N=L/a.
%   The output from c=SYMPHASE(dgt(f,g,a,M),a) is given by
%
%                  L-1 
%     c(m+1,n+1) = sum f(l+1)*exp(-2*pi*i*m*(l-n*a/2)/M)*conj(g(l-a*n+1)), 
%                  l=0  
%
%   where m=0,...,M-1 and n=0,...,N-1 and l-an is computed modulo L.
%
%   SYMPHASE(c,a,'lt',lt) does the same for a non-separable lattice
%   specified by lt. Please see the help of MATRIX2LATTICETYPE for a
%   precise description of the parameter lt.
%
%   See also: dgt, phaselock, phaseunlock
%
%   References:
%     E. Chassande-Mottin, I. Daubechies, F. Auger, and P. Flandrin.
%     Differential reassignment. Signal Processing Letters, IEEE,
%     4(10):293-294, 1997.
%     
%
%   Url: http://ltfat.github.io/doc/gabor/symphase.html

% Copyright (C) 2005-2015 Peter L. Soendergaard <peter@sonderport.dk>.
% This file is part of LTFAT version 2.1.1
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   AUTHORS : Peter Balazs, Peter L. Søndergaard.

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.lt=[0 1];
[flags,kv]=ltfatarghelper({},definput,varargin);

if  (prod(size(a))~=1 || ~isnumeric(a))
  error('a must be a scalar');
end;

if rem(a,1)~=0
  error('a must be an integer');
end;

M=size(c,1);
N=size(c,2);
L=N*a;
b=L/M;

if rem(b,1)~=0
  error('Lattice error. The a parameter is probably incorrect.');
end;

TimeInd = (0:(N-1))*a;
FreqInd = (0:(M-1));

phase = FreqInd'*TimeInd;
phase = mod(phase,M);
phase = exp(1i*pi*phase/M);

if kv.lt(1)>0 
    % truly non-separable case
    
    for n=0:(N-1)
        w = mod(n*kv.lt(1)/kv.lt(2),1);
        phase(:,n+1) = phase(:,n+1)*exp(pi*1i*a*w*n/M);
    end
end

% Handle multisignals
c=bsxfun(@times,c,phase);


