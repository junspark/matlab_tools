function [w] = mls3dweight(r,rho_I)

if r>rho_I
   w=0;
else    
   a=r/rho_I;
   w=1-6*a^2+8*a^3-3*a^4;
end

