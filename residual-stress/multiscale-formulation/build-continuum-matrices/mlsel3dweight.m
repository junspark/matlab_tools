function [w] = mlsel3dweight(vec,orient,a,b,c)
%
% Weight is determined from the ellipse equation
% (x/a)^2+(y/b)^2+(z/c)^2=1 on the ellipsoid
% (x/a)^2+(y/b)^2+(z/c)^2>1 outside the ellipsoid
%
% View the vector in the reference of ellipsoid
xyz=orient'*vec;
% 'vecel' contains the coordinates of the points in the reference of the
% ellipsoid
%
r=sqrt((xyz(1)/a)^2+(xyz(2)/b)^2+(xyz(3)/c)^2);
d=r;
if r>1
   w=0;
%
% Quadratic spline 
% elseif r<0.5
%    w=2/3-(4*d^2)+(4*d^3);
% else
%    w=4/3-(4*d)+(4*d^2)-(4*d^3/3);
%
% Quartic spline
else
   w=1-6*d^2+8*d^3-3*d^4;
% Exponential
% else
%alpha=0.4;
%    w=exp(-d/alpha)
end
%
%
%