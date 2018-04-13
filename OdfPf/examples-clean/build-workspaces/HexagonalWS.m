%
%  Generate workspaces for fundamental region.
%
%-------------------- User Input
%
wsname = 'wshex20x';
%
frbase = HexBaseMesh;
sym    = HexSymmetries;
%
fr_refine   = 20;       % refinement level on FR
sp_refine   = 20;        % refinement level on sphere
per_fiber   = 1000;     % points per fiber
%
% pf_hkls ={[1 0 0], [0 0 1], [1 0 1]};
pf_hkls ={};
%
wsopts = {...
    'MakePoleFigures',   pf_hkls, ...
    'PointsPerFiber',  per_fiber, ...
    'MakeFRL2IP',           'on', ...
    'MakeFRH1IP',           'on', ...
    'MakeSphL2IP',         'off', ...
    'MakeSphH1IP',         'off'  ...
	 };
%
%-------------------- Build workspace
%
fr  = RefineMesh(frbase, fr_refine, sym);
fr.symmetries = sym;
%
sph0    = SphBaseMesh(2, 'Hemisphere', 'off'); 
sph     = RefineMesh(sph0, sp_refine);
sph.crd = UnitVector(sph.crd);
%
wshex   = Workspace(fr, sph, wsopts{:});

%%%
% nHarm       = 109;      % number of harmonics
% [dh, eVals] = DiscreteHarmonics(wshex.frmesh, nHarm);
% wshex.frmesh.dh     = dh;
% wshex.frmesh.eVals  = eVals;
% plot(eVals)
% % return
% save wshex12x wshex wsopts

PlotFR(wshex.frmesh, ones(wshex.frmesh.numind,1), 'ShowMesh', 'on', 'Symmetries', 'Hexagonal')