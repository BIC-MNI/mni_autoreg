% return the value of the 2nd derivative of the data, based on cubic
% interpolation through the four equally spaced data function value
% points.
% xpos is the fractional distance from the x position of f2 for which the
% value is to be returned.
function v = cubic_pp(f1,f2,f3,f4,xpos)

%  for f1 = f(-1),
%      f2 = f(0),
%      f3 = f(1),
%      f4 = f(2), 
%
%      a1 =                f2;
%      a2 =  -1.0*f1/3.0 - f2/2.0 + f3     - f4/6.0;
%      a3 =       f1/2.0 - f2     + f3/2.0 ;
%      a4 =  -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
%
%  to estimate function f() at xpos, using cubic interpolation:
%
%      f   = a1 + a2*xpos + a3*xpos^2 + a4*xpos^3
%      fp  = a2 + 2*a3*xpos + 3*a4*xpos*xpos;      first derivitive
%      fpp = 2*a3 + 6*a4*xpos;	                   second derivitive 

a3 =      f1/2.0 - f2     + f3/2.0 ;
a4 = -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;

v =  2*a3 + 6*a4*xpos;