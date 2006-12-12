% Example program of the optimal interpolation toolbox


% number of observations to interpolate

on = 200;

% create randomly located observations within 
% the square [0 1] x [0 1]

x = rand(1,on);
y = rand(1,on);

% the underlying function to interpolate

f = sin(6*x) .* cos(6*y);

% the error variance of the observations

var = 0.0 * ones(on,1);

% the grid onto which the observations are interpolated

[xi,yi] = ndgrid(linspace(0,1,100));

% the correlation length in x and y direction

lenx = 0.1;
leny = 0.1;

% number of influential observations

m = 30;

% run the optimal interpolation
% fi is the interpolated field and vari is its error variance

[fi,vari] = optiminterp2(x,y,f,var,lenx,leny,m,xi,yi);

