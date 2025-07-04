function filter = gfb_filter_new(arg1,arg2,arg3,arg4,arg5)
%GFB_FILTER_NEW   Constructor of a cascaded gammatonefilter
%   Usage:  gfb_filter_new(a_tilde,gamma_filter_order);
%           gfb_filter_new(fs,fc,gamma_filter_order,bandwidth_factor);
%           gfb_filter_new(fs,fc,bw,attenuation_db,gamma_filter_order);
%
%   Input arguments:
%      a_tilde             : Complex filter constant
%      gamma_filter_order  : Gammatone filter order
%      fs                  : Sampling rate
%      fc                  : Centre frequency
%      a_tilde             : Complex filter constant
%      bandwidth_factor    : Bandwidth factor, default value is 1.0 
%      bw                  : Filter bandwidth (in Hz)
% 
%   GFB_FILTER_NEW(a_tilde,gamma_filter_order) specifies the complex
%   filter coefficients directly.
%
%   GFB_FILTER_NEW(fs,fc,gamma_filter_order,bandwidth_factor) computes
%   filter coefficient from sampling rate, center frequency, and order of
%   the gammatone filter.  The filters will have 1 ERB equivalent
%   rectangular bandwidth, times bandwidth_factor.  Filter coefficient are
%   computed from equations (13),(14)[Hohmann 2002].
%
%   GFB_FILTER_NEW(fs,fc,bw,attenuation_db,gamma_filter_order) Computes
%   the filter coefficients from the sampling rate, center frequency, the
%   desired bandwidth with respect to the given attenuation, and the order
%   of the gammatone filter. Filter coefficient are computed as in equations
%   (11),(12)[Hohmann 2002] (section 2.3).
%
%   References:
%     V. Hohmann. Frequency analysis and synthesis using a gammatone
%     filterbank. Acta Acustica united with Acoustica, 88(3):433-442, 2002.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/doc/filters/gfb_filter_new.php

% Copyright (C) 2009-2015 Piotr Majdak and Peter L. Søndergaard.
% This file is part of AMToolbox version 0.9.7
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

%   AUTHOR: tp


filter.type = 'gfb_Filter';
if (nargin == 2)
  filter.coefficient = arg1;
  filter.gamma_order = arg2;
elseif (nargin == 3) || (nargin == 4)
  fs    = arg1;
  fc = arg2;
  filter.gamma_order  = arg3;
  bandwidth_factor    = 1.0;
  if (nargin == 4)
    bandwidth_factor  = arg4;
  end
  global GFB_L;
  global GFB_Q;
  gfb_set_constants;

  % equation (13) [Hohmann 2002]:
  audiological_erb = (GFB_L + fc / GFB_Q) * bandwidth_factor;
  % equation (14), line 3 [Hohmann 2002]:
  a_gamma          = (pi * factorial(2*filter.gamma_order - 2) * ...
                      2 ^ -(2*filter.gamma_order - 2) /              ...
                      factorial(filter.gamma_order - 1) ^ 2);
  % equation (14), line 2 [Hohmann 2002]:
  b                = audiological_erb / a_gamma;
  % equation (14), line 1 [Hohmann 2002]:
  lambda           = exp(-2 * pi * b / fs);
  % equation (10) [Hohmann 2002]:
  beta             = 2 * pi * fc / fs;
  % equation (1), line 2 [Hohmann 2002]:
  filter.coefficient   = lambda * exp(1i * beta);
elseif (nargin == 5)
  fs    = arg1;
  fc = arg2;
  bw        = arg3;
  attenuation_db      = arg4;
  filter.gamma_order  = arg5;

  % equation (12), line 4 [Hohmann 2002]:
  phi    =  pi * bw / fs;
  % equation (12), line 3 [Hohmann 2002]:
  u      = -attenuation_db/filter.gamma_order;
  % equation (12), line 2 [Hohmann 2002]:
  p      =  (-2 + 2 * 10^(u/10) * cos(phi)) / (1 - 10^(u/10));
  % equation (12), line 1 [Hohmann 2002]:
  lambda = -p/2 - sqrt(p*p/4 - 1);
  % equation (10) [Hohmann 2002]:
  beta   =  2*pi*fc/fs;
  % equation (1), line 2 [Hohmann 2002]:
  filter.coefficient   = lambda * exp(1i*beta);
else
  error ('gfb_filter_new needs either 2, 3, 4 or 5 arguments');
end

% normalization factor from section 2.2 (text) [Hohmann 2002]:
filter.normalization_factor = ...
    2 * (1 - abs(filter.coefficient)) ^ filter.gamma_order;

filter.state = zeros(1, filter.gamma_order);

