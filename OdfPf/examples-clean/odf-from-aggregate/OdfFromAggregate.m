%
%  ODF From Aggregate.
%
%  This script illustrates the methods available for
%  constructing an ODF from a collection of individual 
%  points.
%
%
%-------------------- User Input
%
DegToRad = (pi/180);     % matlab `deg2rad' requires toolbox
wsname = 'wscub';        % workspace name
numpts = 2000;           % number of points
std_agg = DegToRad * 5;  % std deviation for AggregateFunction
std_sm  = DegToRad * 10; % std deviation for smoothing
%
%-------------------- Execution
%
%  Load workspace for fundamental region.
%
addpath('../build-workspaces');
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname);
%
%  Generate a random aggregate.
%
wts  = ones(1, numpts);
quat = UnitVector(randn(4, numpts));
rod  = ToFundamentalRegion(quat, ws.frmesh.symmetries);
%
%  Method 1:  DiscreteDelta  (fast)
%
tic
[elem, ecrd] = MeshCoordinates(ws.frmesh, rod);
odf1 = DiscreteDelta(ws.frmesh, ws.frmesh.l2ip, elem, ecrd, wts);
odf1 = odf1./MeanValue(odf1, ws.frmesh.l2ip);
t = toc;
disp(['Time for DiscreteDelta:  ', num2str(t)]);
%
%  Method 1a:  Smoothed DiscreteDelta (slow)
%
gqrule = QRuleGlobal(ws.frmesh, ws.frmesh.qrule, @RodMetric);
tic
odf1a = SmoothFunction(odf1, ws.frmesh, ...
		       ws.frmesh.qrule.pts, gqrule, ...
		       @RodGaussian, std_sm, ws.frmesh.symmetries);
odf1a = odf1a ./MeanValue(odf1a, ws.frmesh.l2ip);
t = toc;
disp(['Time for SmoothFunction:  ', num2str(t)]);
%
%  Method 2:  AggregateFunction (medium)
%
pts   = ws.frmesh.crd(:, 1:ws.frmesh.numind);
agg   = rod;
wts   = ones(1, numpts);
%
PointFun = @RodGaussian;
sym      = CubSymmetries;
%
tic
odf2 = AggregateFunction(pts, agg, wts, PointFun, std_agg, sym);
odf2 = odf2./MeanValue(odf2, ws.frmesh.l2ip);
t = toc;
disp(['Time for AggregateFunction:  ', num2str(t)]);
%
%  Create DX output files.
%
Ndata = {'odf1', odf1', 'odf1a', odf1a', 'odf2', odf2'};
ExportDX('agg-odf', ws.frmesh, Ndata);
%
save odfs
