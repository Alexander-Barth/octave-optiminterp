%% Copyright (C) 2006 Alexander Barth <barth.alexander@gmail.com>
%%
%% This program is free software; you can redistribute it and/or modify it under
%% the terms of the GNU General Public License as published by the Free Software
%% Foundation; either version 3 of the License, or (at your option) any later
%% version.
%%
%% This program is distributed in the hope that it will be useful, but WITHOUT
%% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
%% details.
%%
%% You should have received a copy of the GNU General Public License along with
%% this program; if not, see <http://www.gnu.org/licenses/>.

% Tests 1D, 2D, 3D and 4D optimal interpolation.
% All tests should pass; any error indicates that either 
% there is a bug in the optimal interpolation package or 
% that it is incrorrectly installed.



printf('Benchmark 3D-optimal interpolation: ');


% grid of background field
ni = 50
[xi,yi,zi] = ndgrid(linspace(0,1,ni));
fi_ref = cos(6*xi) .* cos(6*yi) .* cos(6*zi);

% grid of observations
[x,y,z] = ndgrid(linspace(0,1,10));
x = x(:);
y = y(:);
z = z(:);

on = numel(x);
var = 0.01 * ones(on,1);
f = cos(6*x) .* cos(6*y) .* cos(6*z);

m = 20;

tic
  [fi,vari] = optiminterp3(x,y,z,f,var,0.1,0.1,0.1,m,xi,yi,zi);
  ctime = toc


  % ctime =  12.863