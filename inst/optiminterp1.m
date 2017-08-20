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
## @deftypefn {Loadable Function} {[@var{fi},@var{vari}]} = optiminterp1(@var{x},@var{f},@var{var},@var{lenx},@var{m},@var{xi})
## Performs a local 1D-optimal interpolation (objective analysis).
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

  [fi,vari] = optiminterpn(x,f,var,lenx,m,xi);

endfunction

%!test
%! # grid of background field
%! xi = linspace (0, 1, 50)';
%! fi_ref(:,1) = sin (xi * 6);
%! fi_ref(:,2) = cos (xi * 6);
%!
%! # grid of observations
%! x = linspace (0, 1, 20)';
%!
%! on = numel (x);
%! var = 0.01 * ones (on, 1);
%! f(:,1) = sin (x * 6);
%! f(:,2) = cos (x * 6);
%!
%! m = 15;
%!
%! [fi, vari] = optiminterp1 (x, f, var, 0.1, m, xi);
%!
%! rms = sqrt (mean ((fi_ref(:) - fi(:)) .^ 2));
%!
%! assert (rms <= 0.005, "unexpected large difference with reference field");
