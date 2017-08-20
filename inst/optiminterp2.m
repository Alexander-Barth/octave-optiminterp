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
## @deftypefn {Loadable Function} {[@var{fi},@var{vari}] = } optiminterp2(@var{x},@var{y},@var{f},@var{var},@var{lenx},@var{leny},@var{m},@var{xi},@var{yi})
## Performs a local 2D-optimal interpolation (objective analysis).
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

function [fi,vari] = optiminterp2(x,y,f,var,lenx,leny,m,xi,yi)

  [fi,vari] = optiminterpn(x,y,f,var,lenx,leny,m,xi,yi);

endfunction

%!demo
%! # the grid onto which the observations are interpolated
%! [xi, yi] = ndgrid (linspace (0, 1, 100));
%!
%! # background estimate or first guess
%! xb = 10 + xi;
%!
%! # number of observations to interpolate
%! on = 200;
%!
%! # create randomly located observations within 
%! # the square [0 1] x [0 1]
%! x = rand (1, on);
%! y = rand (1, on);
%!
%! # the underlying function to interpolate
%! yo = 10 + x + sin (6 * x) .* cos (6 * y);
%!
%! # the error variance of the observations divided by the error 
%! # variance of the background field
%! var = 0.1 * ones (on, 1);
%!
%! # the correlation length in x and y direction
%! lenx = 0.1;
%! leny = 0.1;
%!
%! # number of influential observations
%! m = 30;
%!
%! # subtract the first guess from the observations
%! # (DON'T FORGET THIS - THIS IS VERY IMPORTANT)
%! Hxb = interp2 (xi(:,1), yi(1,:), xb', x, y);
%! f = yo - Hxb;
%!
%! # run the optimal interpolation
%! # fi is the interpolated field and vari is its error variance
%! [fi, vari] = optiminterp2 (x, y, f, var, lenx, leny, m, xi, yi);
%!
%! # Add the first guess back
%! xa = fi + xb;
