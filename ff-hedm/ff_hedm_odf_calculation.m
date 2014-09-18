clear all
close all
clc

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
wsname = 'wscub4x';        % workspace name
numpts = 2000;           % number of points
std_agg = DegToRad * 5;  % std deviation for AggregateFunction
std_sm  = DegToRad * 10; % std deviation for smoothing
%
%-------------------- Execution
%
%  Load workspace for fundamental region.
%
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname)

% grains	= load('.\hedm_HTUPS_unirr_s1_MultiRing\Layer1_ring1_t70_2_t50\Grains.csv');    %%%%%%%% PROMISING
grains	= load('.\hedm_HTUPS_irr_s1_MultiRing\Layer1\Grains.csv');          %%%%%%%% PROMISING

numpts  =  size(grains,1);
wts     = ones(1, numpts);
%
%  Generate a random aggregate.
%
RESRF2APS   = RMatOfQuat(QuatOfESRF2APS);
% RESRF2APS   = eye(3,3);

nGrains     = size(grains, 1);
for i = 1:1:nGrains
    RMats(:,:,i)   =  RESRF2APS*reshape(grains(i,2:10), 3, 3)';
end
qsym    = CubSymmetries; Rsym    = RMatOfQuat(qsym);
quat    = QuatOfRMat(RMats);
rod     = ToFundamentalRegion(quat, ws.frmesh.symmetries);

tic
[elem, ecrd] = MeshCoordinates(ws.frmesh, rod);
odf1 = DiscreteDelta(ws.frmesh, ws.frmesh.l2ip, elem, ecrd, wts);
odf1 = odf1./MeanValue(odf1, ws.frmesh.l2ip);
t = toc;
disp(['Time for DiscreteDelta:  ', num2str(t)]);
%
%  Method 1a:  Smoothed DiscreteDelta (slow)
%
% gqrule = QRuleGlobal(ws.frmesh, ws.frmesh.qrule, @RodMetric);
% tic
% odf1a = SmoothFunction(odf1, ws.frmesh, ...
% 		       ws.frmesh.qrule.pts, gqrule, ...
% 		       @RodGaussian, std_sm, ws.frmesh.symmetries);
% odf1a = odf1a ./MeanValue(odf1a, ws.frmesh.l2ip);
% t = toc;
% disp(['Time for SmoothFunction:  ', num2str(t)]);
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
% Ndata = {'odf1', odf1', 'odf1a', odf1a', 'odf2', odf2'};
% ExportDX('agg-odf', ws.frmesh, Ndata);
%
% save odfs
