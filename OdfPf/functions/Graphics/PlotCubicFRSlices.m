function PlotCubicFRSlices(mesh, odf)
% PLOTCUBICFRSLICES - Plot slices of ODF on coordinate planes.
%   
%   PlotCubicFRSlices(mesh, odf)
%
%   mesh is a MeshStructure on the fundamental region
%
%   odf is a vector, the function values to plot on the 
%                    reduced (nonequivalent) set of nodes
% 
crd = mesh.crd;  % Unpack input structures.
con = mesh.con;
eqv = mesh.eqv;
%
odf = ToAllNodes(odf, eqv);
%
%-------------------- Generate and plot mesh for slices.
%
xmax  = 0.4142; xmin = -xmax;
ndiv  = 10; npts = (2*ndiv+1)^2;
range = xmax/ndiv;
%
[x, y] = meshgrid(range*(-ndiv:ndiv));
x = reshape(x, 1, npts); 
y = reshape(y, 1, npts);
z = zeros(1, npts);
xycon = delaunay(x, y); % numtri x 3
%
slice_crd = [...
    x y z; ...
    y z x; ...
    z x y; ...
    ];
slice_con = [...
    xycon; xycon + npts; xycon + 2*npts; ...
    ]';
slice_mesh = MeshStructure(slice_crd, slice_con, []);
%
%  Find slice locations in 3D mesh.
%
[e, ecrd] = tsearchn(crd', con', slice_crd');
%fails = find(~finite(e));
%%
%fail_crds = slice_crd(fails,:)
%%
eodf      = odf(con(:, e))';
slice_odf = dot(ecrd, eodf, 2);
%
PlotSurface(slice_mesh, slice_odf);
%
%-------------------- Show edges.
%
zedge = [xmin xmin 0; ...
	 xmax xmin 0; ...
	 xmax xmax 0; ...
	 xmin xmax 0; ...
	 xmin xmin 0];
%
linex = [zedge(:,1) zedge(:,2) zedge(:,3)];
liney = [zedge(:,2) zedge(:,3) zedge(:,1)];
linez = [zedge(:,3) zedge(:,1) zedge(:,2)];
l = line(linex, liney, linez, 'Color', 'k');
%
%  Show coordinate axes.
%
axesx = [...
    xmin  0     0 ; ...
    xmax  0     0];
axesy = [...
    0     xmin  0 ; ...
    0     xmax  0];
axesz = [...
    0     0     xmin ; ...
    0     0     xmax];
%
l = line(axesx, axesy, axesz, 'Color', 'k');
%
axis off
