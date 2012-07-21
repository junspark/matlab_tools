%
%  Spherical harmonics.
%
dim = 3;  %  dimension of sphere
div = 5;  %  number of subdivisions
nev = 20; %  number of eigenvalues
MeshType = 'hfr'; % 'sphere', 'cfr', 'hfr'
%
%  Quadrature rule.
%
if (dim == 2)
  qrname = 'qr_trid05p07';  
else
  qrname = 'qr_tetd05p14';
end
%
%  Base mesh.
%
if strcmp(MeshType, 'sphere')
  %
  %  :  S^2,  S^3
  %
  Snbase = SphBaseMesh(dim);
  Snmesh = RefineMesh(Snbase, div);
  Snmesh.crd = UnitVector(Snmesh.crd);
elseif strcmp(MeshType, 'cfr')
  %
  %  :  CFR
  %
  Cbase = CubBaseMesh();
  Cmesh  = RefineMesh(Cbase, div, CubSymmetries);
  Snmesh = Cmesh;
  Snmesh.crd = QuatOfRod(Snmesh.crd);
  %
elseif strcmp(MeshType, 'hfr')
  %
  %  :  HFR
  %
  Hbase  = HexBaseMesh;
  Hmesh  = RefineMesh(Hbase, div, HexSymmetries);
  Snmesh = Hmesh;
  Snmesh.crd = QuatOfRod(Snmesh.crd);
  %
end
%
%
%--------------------*--------------------------------------------------
%
load(qrname);
eval(['qr = ', qrname, ';']);
%
disp('evaluating matrices ...'); tic
h1sip = SphH1SIP(Snmesh, qr);
l2ip  = SphL2IP(Snmesh, qr);
disp('done'); toc
%
disp('computing eigenvalues ...'); tic
OPTS.disp = 0;
[V, D] = eigs(h1sip, l2ip, nev, 'SM', OPTS);
d = sort(diag(D))';
V = V(:, nev:-1:1);
disp('done'); toc
