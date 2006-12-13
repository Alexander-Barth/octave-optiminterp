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
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Loadable Function} {[@var{fi},@var{vari}] = } optiminterp(@var{x},@var{y},@var{f},@var{var},@var{lenx},@var{leny},@var{m},@var{xi},@var{yi})
## Performs a 2D-optimal interpolation (objective analysis).
##
## Every elements in @var{f} corresponds to a data point (observation)
## at location @var{x},@var{y} with the error variance @var{var}.
##
## @var{lenx} and @var{leny} are correlation length in x-direction
## and y-direction respectively. 
## @var{m} represents the number of influential points.
##
## @var{xi} and @var{yi} are the data points where the field is
## interpolated. @var{fi} is the interpolated field and @var{vari} is 
## its error variance.
##
## The background field of the optimal interpolation is zero.
## For a different background field, the background field
## must be subtracted from the observation, the difference 
## is mapped by OI onto the background grid and finally the
## background is added back to the interpolated field.
## The error variance of the background field is assumed to 
## have a error variance of one.
## @end deftypefn

## Copyright (C) 2006, Alexander Barth
## Author: Alexander Barth <abarth@marine.usf.edu>

function [fi,vari] = optiminterp2(x,y,f,var,lenx,leny,m,xi,yi)

if (isscalar(var))
  var = var*ones(size(x));
end

if isvector(f) & size(f,1) == 1
   f = f';
end

on = numel(x);
nf = size(f,3);

if (on*nf ~= numel(f) & on ~= numel(y)  & on ~= numel(var))
  error('optiminterp2: x,y,f,var must have the same number of elements');
  return
end

if (numel(xi) ~= numel(yi))
  error('optiminterp2: xi and yi must have the same number of elements');
end

gx(:,1) = xi(:); 
gx(:,2) = yi(:); 

ox(:,1) = x(:);
ox(:,2) = y(:);

f=reshape(f,[on nf]);

%whos ox gx f

[fi,vari] = optiminterp(ox,f,var,[lenx leny],m,gx);

fi = reshape(fi,[size(xi) nf]);
vari = reshape(vari,size(xi));
