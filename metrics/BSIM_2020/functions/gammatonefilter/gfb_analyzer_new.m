function analyzer = gfb_analyzer_new(fs,flow,basef,fhigh,filters_per_ERBaud,gamma_order,bandwidth_factor)
%GFB_ANALYZER_NEW  Construct new analyzer object
%   Usage:  analyzer = gfb_analyzer_new(fs,flow, basef, fhigh,filters_per_ERBaud,gamma_order,bandwidth_factor)
%
%   Input parameters:
%      fs                 : The sampling frequency of the signals on which
%                           the analyzer will operate
%      flow               : The lowest possible center frequency of a
%                           contained gammatone filter
%      basef              : "base frequency". One of the gammatone filters of the analyzer
%                           will have this center frequency.  Must be >= flow
%      fhigh              : The highest possible center frequency of a
%                           contained gammatone filter.  Must be >=  basef
%      filters_per_ERBaud : The density of gammatone filters on the ERB
%                           scale.
%      gamma_order        : The order of the gammatone filters in this
%                           filterbank. If unspecified, the default value from
%                           gfb_set_constants.m is used.
%      bandwidth_factor   : The bandwidth parameter of the individual filters
%                           is calculated from the Equivalent Rectangular
%                           Bandwidth (ERB) according to equation 14 in
%                           Hohmann (2002). ERB is taken from the Glasberg &
%                           Moore formula for a specific center frequency
%                           (equation 13 in Hohmann (2002)).
%                           Using this parameter, it is possible to widen or
%                           narrow all filters of the filterbank with a
%                           constant bandwidth factor.
%                           Default value is 1.0
%
%   Output parameters:
%      analyzer           : The constructed gfb_analyzer object.
%    
%   GFB_ANALYZER_NEW constructs a new gfb_analyzer object.  The analyzer
%   implements the analysis part of a gammatone filterbank as described
%   in Hohmann (2002).
%
%   It consists of several all-pole gammatone filters; each
%   one with a bandwidth of 1 ERBaud (times bandwidth_factor),
%   and an order of gamma_order.
%
%   The center frequencies of the individual filters are computed as
%   described in section 3 of Hohmann (2002).
%
%
%   Url: http://amtoolbox.sourceforge.net/doc/filters/gfb_analyzer_new.php

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
  
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Jan, Sep 2003, Nov 2006, Jan 2007

if (nargin < 6)
  % The order of the gammatone filter is derived from the global constant
  % GFB_PREFERED_GAMMA_ORDER defined in "gfb_set_constants.m".  Usually,
  % this is equal to 4.
  global GFB_PREFERED_GAMMA_ORDER;
  gfb_set_constants;
  gamma_order = GFB_PREFERED_GAMMA_ORDER;
end
if (nargin < 7)
  bandwidth_factor = 1.0;
end

% To avoid storing information in global variables, we use Matlab
% structures:
analyzer.type                          = 'gfb_analyzer';
analyzer.fs         = fs;
analyzer.flow     = flow;
analyzer.basef = basef;
analyzer.fhigh     = fhigh;
analyzer.filters_per_ERBaud            = filters_per_ERBaud;
analyzer.bandwidth_factor              = bandwidth_factor;
analyzer.fast                          = 0;


%
analyzer.center_frequencies_hz = ...
    erbspacebw(flow,fhigh,1/filters_per_ERBaud,basef);

% This loop actually creates the gammatone filters:
for band = [1:length(analyzer.center_frequencies_hz)]
  center_frequency_hz = analyzer.center_frequencies_hz(band);

  % Construct gammatone filter with one ERBaud bandwidth:
  analyzer.filters(1,band) = ...
      gfb_filter_new(fs, center_frequency_hz, ...
                     gamma_order, bandwidth_factor);
end

