function y = pfunc_ka12(p,x)
bk(1)   = p(1);
bk(2)   = p(2);

pg1     = p(3:5);
pg2     = p(6:8);
pl1     = [p(9) p(10) p(5)];
pl2     = [p(11) p(12) p(8)];
n1      = p(13);
n2      = p(14);

ybk = fcnPolynomial(bk,x);

g1  = pfunc_Gaussian(pg1,x);
l1  = pfunc_Lorentzian(pl1,x);

g2  = pfunc_Gaussian(pg2,x);
l2  = pfunc_Lorentzian(pl2,x);

y1  = n1*g1 + (1-n1)*l1;
y2  = n2*g2 + (1-n2)*l2;

y   = ybk + y1 + y2;