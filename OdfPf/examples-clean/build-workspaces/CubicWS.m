%
%  Generate workspaces for cubic fundamental region.
%
%-------------------- User Input
%
fr_refine   =    4;  % refinement level on FR
sp_refine   =    15;  % refinement level on sphere
per_fiber   = 1000;  % points per fiber
nHarm       = 130;   % number of harmonics
%
pf_hkls ={[1 1 1], [2 0 0], [2 2 0]};
%
wsopts = {...
    'MakePoleFigures',   pf_hkls, ...
    'PointsPerFiber',  per_fiber, ...
    'MakeFRL2IP',           'on', ...
    'MakeFRH1IP',           'on', ...
    'MakeSphL2IP',          'on', ...
    'MakeSphH1IP',          'on'  ...
	 };

%
%-------------------- Build workspace
%
cfr0 = CubBaseMesh;
csym = CubSymmetries;
cfr  = RefineMesh(cfr0, fr_refine, csym);
cfr.symmetries = csym;
%
sph0 = SphBaseMesh(2, 'Hemisphere', 'off'); % 2d sphere
sph  = RefineMesh(sph0, sp_refine);
sph.crd = UnitVector(sph.crd);

%
wscub = Workspace(cfr, sph, wsopts{:});

[dh, eVals] = DiscreteHarmonics(wscub.frmesh, nHarm);
wscub.frmesh.dh     = dh;
wscub.frmesh.eVals  = eVals;

save wscub wscub wsopts

PlotFR(wscub.frmesh, ones(wscub.frmesh.numind,1), 'ShowMesh', 'on')