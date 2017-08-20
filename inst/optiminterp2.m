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

%!xtest
%! # Tests 2D optimal interpolation and compares result 
%! # (analysis and error field) to global OI solution
%!
%! # grid of background field
%! [xi, yi] = ndgrid (linspace (0, 1, 30));
%! fi_ref = sin (xi * 6) .* cos(yi * 6);
%!
%! # grid of observations
%! [x, y] = ndgrid (linspace (0, 1, 6));
%! x = x(:);
%! y = y(:);
%!
%! on = numel (x);
%! var = 0.01 * ones (on, 1);
%! f = sin (x * 6) .* cos (y * 6);
%!
%! len = 0.1;
%! m = min (30, on);
%!
%! # covariance function
%! # gaussian
%! #bcovar2 = @(d2) exp (-d2 / len ^ 2) ;
%! # diva
%! bcovar2 = @(d2) max (sqrt (d2) / len, eps) .* besselk (1, max (sqrt (d2) / len, eps));
%!
%! # P: covariance between grid points (xi,yi) and grid points (xi,yi)
%! P = zeros (numel (xi), numel (xi));
%!
%! for j = 1:numel (xi)
%!   for i = 1:numel (xi) 
%!     P(i,j) = (xi(i) - xi(j)) ^ 2 + (yi(i) - yi(j)) ^ 2;
%!   endfor
%! endfor
%! P = bcovar2 (P);
%!
%! # HPH: covariance between observation points (x,y) and observation points (x,y)
%! HPH = zeros (numel (x), numel (x));
%!
%! for j = 1:numel (x)
%!   for i = 1:numel (x)
%!     HPH(i,j) = (x(i) - x(j)) ^ 2 + (y(i) - y(j)) ^ 2;
%!   endfor
%! endfor
%! HPH = bcovar2 (HPH);
%!
%! # PH: covariance between grid points (xi,yi) and observation points (x,y)
%! PH = zeros (numel (xi),numel (x));
%!
%! for j = 1:numel (x)
%!   for i = 1:numel (xi) 
%!     PH(i,j) = (xi(i) - x(j)) ^ 2 + (yi(i) - y(j)) ^ 2;
%!   endfor
%! endfor
%! PH = bcovar2 (PH);
%!
%! R = diag (var);
%!
%! # call optiminterp
%! [fi, vari] = optiminterp2 (x, y, f, var, len, len, m, xi, yi);
%!
%! # Kalman gain
%! K = PH * inv (HPH + R);
%!
%! # analysis
%! fi2 = K * f;
%!
%! # error field
%! vari2 = diag (P - K * PH');
%!
%! # transform vectors into 2d-arrays
%! fi2 = reshape (fi2, size (fi));
%! vari2 = reshape (vari2, size (fi));
%!
%! rms = sqrt (mean ((fi2(:) - fi(:)) .^ 2));
%!
%! assert (rms <= 1e-4, "unexpected large RMS difference (analysis)");
%!
%! rms = sqrt (mean ((vari2(:) - vari(:)) .^ 2));
%!
%! assert (rms <= 1e-4, "unexpected large RMS difference (error field)");

%!test
%! # grid of background field
%! [xi, yi] = ndgrid (linspace (0, 1, 30));
%!
%! fi_ref(:,:,1) = sin (xi * 6) .* cos (yi * 6);
%! fi_ref(:,:,2) = cos (xi * 6) .* sin (yi * 6);
%!
%! # grid of observations
%! [x, y] = ndgrid (linspace (0, 1, 20));
%!
%! on = numel (x);
%! var = 0.01 * ones (on, 1);
%! f(:,:,1) = sin (x * 6) .* cos (y * 6);
%! f(:,:,2) = cos (x * 6) .* sin (y * 6);
%!
%! m = 30;
%!
%! [fi, vari] = optiminterp2 (x, y, f, var, 0.1, 0.1, m, xi, yi);
%!
%! rms = sqrt (mean ((fi_ref(:) - fi(:)) .^ 2));
%!
%! assert (rms <= 0.005, "unexpected large difference with reference field");
