## Copyright (C) 2006 Alexander Barth <barth.alexander@gmail.com>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Loadable Function} {[@var{fi},@var{vari}] = } optiminterp3(@var{x},@var{y},@var{z},@var{f},@var{var},@var{lenx},@var{leny},@var{lenz},@var{m},@var{xi},@var{yi},@var{zi})
## Performs a local 3D-optimal interpolation (objective analysis).
##
## Every elements in @var{f} corresponds to a data point (observation)
## at location  @var{x}, @var{y}, @var{z} with the error variance var
##
## @var{lenx},@var{leny} and @var{lenz} are correlation length in x-,y- and z-direction
## respectively. 
## @var{m} represents the number of influential points.
##
## @var{xi},@var{yi} and @var{zi} are the data points where the field is
## interpolated. @var{fi} is the interpolated field and @var{vari} is 
## its error variance.
##
##
## The background field of the optimal interpolation is zero.
## For a different background field, the background field
## must be subtracted from the observation, the difference 
## is mapped by OI onto the background grid and finally the
## background is added back to the interpolated field.
##
## The error variance of the background field is assumed to 
## have a error variance of one. 
## @end deftypefn

function [fi,vari] = optiminterp3(x,y,z,f,var,lenx,leny,lenz,m,xi,yi,zi)

  [fi,vari] = optiminterpn(x,y,z,f,var,lenx,leny,lenz,m,xi,yi,zi);

endfunction

%!test
%! # grid of background field
%! [xi, yi, zi] = ndgrid (linspace (0, 1, 15));
%!
%! fi_ref(:,:,:,1) = sin (6 * xi) .* cos (6 * yi) .* sin (6 * zi);
%! fi_ref(:,:,:,2) = cos (6 * xi) .* sin (6 * yi) .* cos (6 * zi);
%!
%! # grid of observations
%! [x, y, z] = ndgrid (linspace (0, 1, 10));
%!
%! on = numel (x);
%! var = 0.01 * ones (on, 1);
%! f(:,:,:,1) = sin (6 * x) .* cos (6 * y) .* sin (6 * z);
%! f(:,:,:,2) = cos (6 * x) .* sin (6 * y) .* cos (6 * z);
%!
%! m = 20;
%!
%! [fi, vari] = optiminterp3 (x, y, z, f, var, 0.1, 0.1, 0.1, m, xi, yi, zi);
%!
%! rms = sqrt (mean ((fi_ref(:) - fi(:)) .^ 2));
%!
%! assert (rms <= 0.04, "unexpected large difference with reference field");

%!test
%! # grid of background field
%! [xi, yi, zi] = ndgrid (linspace (0, 1, 15));
%! fi_ref = sin (6 * xi) .* cos (6 * yi) .* sin (6 * zi);
%!
%! # grid of observations
%! [x, y, z] = ndgrid (linspace (0, 1, 10));
%! x = x(:);
%! y = y(:);
%! z = z(:);
%!
%! on = numel (x);
%! var = 0.01 * ones (on, 1);
%! f = sin (6 * x) .* cos (6 * y) .* sin (6 * z);
%!
%! m = 20;
%!
%! [fi, vari] = optiminterp3 (x, y, z, f, var, 0.1, 0.1, 0.1, m, xi, yi, zi);
%!
%! rms = sqrt (mean ((fi_ref(:) - fi(:)) .^ 2));
%!
%! assert (rms <= 0.04, "'unexpected large difference with reference field");
