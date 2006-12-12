% Tests 1D, 2D and 3D optimal interpolation.
% All tests should pass; any error indicates that either 
% there is a bug in the optimal interpolation package or 
% that it is incrorrectly installed.

function test_optiminterp

more off

printf('Testing 1D-optimal interpolation: ');

try
  % grid of background field
  xi = linspace(0,1,50);
  fi_ref = sin( xi*6 );

  % grid of observations
  x = linspace(0,1,20);

  on = numel(x);
  var = 0.01 * ones(on,1);
  f = sin( x*6 );

  m = 15;

  [fi,vari] = optiminterp1(x,f,var,0.1,m,xi);


  rms = sqrt(mean((fi_ref(:) - fi(:)).^2));

  if (rms > 0.005) 
    error('unexpected large difference with reference field');
  end

  disp('OK');

catch
  disp('failed');
  disp(lasterr);
end


printf('Testing 2D-optimal interpolation: ');

try
  % grid of background field
  [xi,yi] = ndgrid(linspace(0,1,30));
  fi_ref = sin( xi*6 ) .* cos( yi*6);

  % grid of observations
  [x,y] = ndgrid(linspace(0,1,20));
  x = x(:);
  y = y(:);

  on = numel(x);
  var = 0.01 * ones(on,1);
  f = sin( x*6 ) .* cos( y*6);

  m = 30;

  [fi,vari] = optiminterp2(x,y,f,var,0.1,0.1,m,xi,yi);


  rms = sqrt(mean((fi_ref(:) - fi(:)).^2));

  if (rms > 0.005) 
    error('unexpected large difference with reference field');
  end

  disp('OK');

catch
  disp('failed');
  disp(lasterr);
end


printf('Testing 3D-optimal interpolation: ');

try
  % grid of background field
  [xi,yi,zi] = ndgrid(linspace(0,1,15));
  fi_ref = sin(6*xi) .* cos(6*yi) .* sin(6*zi);

  % grid of observations
  [x,y,z] = ndgrid(linspace(0,1,10));
  x = x(:);
  y = y(:);
  z = z(:);

  on = numel(x);
  var = 0.01 * ones(on,1);
  f = sin(6*x) .* cos(6*y) .* sin(6*z);

  m = 20;

  [fi,vari] = optiminterp3(x,y,z,f,var,0.1,0.1,0.1,m,xi,yi,zi);


  rms = sqrt(mean((fi_ref(:) - fi(:)).^2));

  if (rms > 0.04) 
    error('unexpected large difference with reference field');
  end

  disp('OK');

catch
  disp('failed');
  disp(lasterr);
end

