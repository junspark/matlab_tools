%
%  Generate workspace for pole figure inversion demo.
%
%-------------------- User Input
%
fr_refine =   3;  % refinement level on FR
per_fiber = 100;  % points per fiber
%
pf_hkls ={[1 1 0], [1 0 0]};
%
pf110 = load('a2mtx2-gage-220.pf');  
pf100 = load('a2mtx2-gage-200.pf');
%
wsopts = {...
    'MakePoleFigures',   pf_hkls, ...
    'PointsPerFiber',  per_fiber, ...
    'MakeFRL2IP',           'on', ...
    'MakeFRH1IP',           'on', ...
    'MakeSphL2IP',          'on', ...
    'MakeSphH1IP',         'off'  ...
	 };
%
%-------------------- Build workspace
%
cfr0 = CubBaseMesh;
csym = CubSymmetries;
cfr  = RefineMesh(cfr0, fr_refine, csym);
cfr.symmetries = csym;
%
ntheta = 72;
nphi   = 18;  % 5 degree increments
sph    = SphCrdMesh(ntheta, nphi);
%
wspfi = Workspace(cfr, sph, wsopts{:});
%
wspfi.pfs(1) = struct('data', pf110);
wspfi.pfs(2) = struct('data', pf100);
%
save wspfi  wspfi wsopts
