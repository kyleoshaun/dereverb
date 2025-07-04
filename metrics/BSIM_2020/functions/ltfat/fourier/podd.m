function f=podd(f,dim)
%PODD   Odd part of periodic function
%   Usage:  fe=podd(f);
%           fe=podd(f,dim);
%
%   PODD(f) returns the odd part of the periodic sequence f.
%
%   PODD(f,dim) does the same along dimension dim.
%
%   See also:  peven, dft, involute, pconv
%
%   Url: http://ltfat.github.io/doc/fourier/podd.html

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
 
if nargin==1 
  f=(f-involute(f))/2;
else
  f=(f-involute(f,dim))/2;
end;

