% Example program of the optimal interpolation toolbox


% the grid onto which the observations are interpolated

[xi,yi] = ndgrid(linspace(0,1,100));

% background estimate or first guess
xb = 10 + xi;

% number of observations to interpolate

on = 200;

% create randomly located observations within 
% the square [0 1] x [0 1]

x = rand(1,on);
y = rand(1,on);

% the underlying function to interpolate

yo = 10 + x + sin(6*x) .* cos(6*y);

% the error variance of the observations divided by the error 
% variance of the background field

var = 0.1 * ones(on,1);

% the correlation length in x and y direction

lenx = 0.1;
leny = 0.1;

% number of influential observations

m = 30;

% subtract the first guess from the observations
% (DON'T FORGET THIS - THIS IS VERY IMPORTANT)

Hxb = interp2(xi(:,1),yi(1,:),xb',x,y);
f = yo - Hxb;

% run the optimal interpolation
% fi is the interpolated field and vari is its error variance

[fi,vari] = optiminterp2(x,y,f,var,lenx,leny,m,xi,yi);

% Add the first guess back

xa = fi + xb;