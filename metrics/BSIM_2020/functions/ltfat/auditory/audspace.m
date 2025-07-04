function [y,bw] = audspace(flow,fhigh,n,varargin)
%AUDSPACE  Equidistantly spaced points on auditory scale
%   Usage: y=audspace(flow,fhigh,n,scale);
%
%   AUDSPACE(flow,fhigh,n,scale) computes a vector of length n*
%   containing values equidistantly scaled on the selected auditory scale
%   between the frequencies flow and fhigh. All frequencies are
%   specified in Hz.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter.
%  
%   [y,bw]=AUDSPACE(...) does the same but outputs the bandwidth between
%   each sample measured on the selected scale.
%  
%   See also: freqtoaud, audspacebw, audfiltbw
%
%   Url: http://ltfat.github.io/doc/auditory/audspace.html

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
  
%   AUTHOR : Peter L. Søndergaard
  
%% ------ Checking of input parameters ---------

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;
  
% Default parameters.

if ~isnumeric(flow) || ~isscalar(flow) 
  error('%s: flow must be a scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) 
  error('%s: fhigh must be a scalar.',upper(mfilename));
end;

if ~isnumeric(n) || ~isscalar(n) || n<=0 || fix(n)~=n
  error('%s: n must be a positive, integer scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;

definput.import={'freqtoaud'};
[flags,kv]=ltfatarghelper({},definput,varargin);


%% ------ Computation --------------------------

audlimits = freqtoaud([flow,fhigh],flags.audscale);

y = audtofreq(linspace(audlimits(1),audlimits(2),n),flags.audscale);

bw=(audlimits(2)-audlimits(1))/(n-1);

% Set the endpoints to be exactly what the user specified, instead of the
% calculated values
y(1)=flow;
y(end)=fhigh;


