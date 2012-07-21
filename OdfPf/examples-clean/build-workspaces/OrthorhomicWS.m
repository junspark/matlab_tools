%
%  Generate workspaces for fundamental region.
%
%-------------------- User Input
%
wsname = 'wsort';
%
frbase = OrtBaseMesh;
sym    = OrtSymmetries;
%
fr_refine =   3;  % refinement level on FR
sp_refine =  10;  % refinement level on sphere
per_fiber = 100;  % points per fiber
%
pf_hkls ={[1 0 0], [0 0 1]};
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
fr  = RefineMesh(frbase, fr_refine, sym);
fr.symmetries = sym;
%
sph0    = SphBaseMesh(2, 'Hemisphere', 'on'); 
sph     = RefineMesh(sph0, sp_refine);
sph.crd = UnitVector(sph.crd);
%
eval([wsname, '= Workspace(fr, sph, wsopts{:});']);
save(wsname, wsname);
