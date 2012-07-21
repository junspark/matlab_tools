%
%  Average properties.
%
mymesh  = LoadMesh('fr-cubic-2');
myqrule = LoadQuadrature('qr_tetd08p43');
gqr     = QRuleGlobal(mymesh, myqrule, @RodMetric);
l2ip    = L2IPMatrix(gqr, NpQpMatrix(mymesh, myqrule));
%
tayfac = load('taylor-factors.dat');
%
wts = sum(l2ip);
avgtayfac = (wts*tayfac)/sum(wts)
