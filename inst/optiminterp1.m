## Copyright (C) 2006 Alexander Barth
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

## -*- texinfo -*-
## @deftypefn {Loadable Function} {[@var{fi},@var{vari}]} = optiminterp(@var{x},@var{f},@var{var},@var{lenx},@var{m},@var{xi})
## Performs a 1D-optimal interpolation (objective analysis).
##
## Every elements in @var{f} corresponds to a data point (observation)
## at location @var{x},@var{y} with the error variance @var{var}.
##
## @var{lenx} is correlation length in x-direction.
## @var{m} represents the number of influential points.
##
## @var{xi} is the data points where the field is
## interpolated. @var{fi} is the interpolated field and @var{vari} is 
## its error variance.
##
## The background field of the optimal interpolation is zero.
## For a different background field, the background field
## must be subtracted from the observation, the difference 
## is mapped by OI onto the background grid and finally the
## background is added back to the interpolated field.
## @end deftypefn



function [fi,vari] = optiminterp1(x,f,var,lenx,m,xi)

if (isscalar(var))
  var = var*ones(size(x));
end

if isvector(f) & size(f,1) == 1
   f = f';
end

if (size(f,1) ~= numel(x)  & numel(f) ~= numel(var))
  error('optiminterp2: x,var must have the same number of elements');
end

%whos ox f var lenx m

gx(:,1) = xi(:);
ox(:,1) = x(:);

%whos ox f var lenx m
[fi,vari] = optiminterp(ox,f,var,lenx,m,gx);


%fi = reshape(fi,[numel(xi) size(f,2)]);
%vari = reshape(vari,size(xi));
