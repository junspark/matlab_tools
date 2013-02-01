function [out] = borders(x,y)
Ro=15;
Ri=6;
out=1;
% check that say the coordinates are outside the borders
if (x>=0) & (y>=0) & (x^2+y^2<=Ro^2) & (x^2+y^2>=Ri^2)
    out=0;
end
