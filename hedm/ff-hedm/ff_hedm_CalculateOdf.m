clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ODF From Aggregate - MIDAS RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DegToRad = (pi/180);     % matlab `deg2rad' requires toolbox
wsname = 'wscub4x';      % workspace name
std_agg = DegToRad * 5;  % std deviation for AggregateFunction
std_sm  = DegToRad * 10; % std deviation for smoothing

pname   = './examples/';
fname   = 'Grains_example.csv';

% ROTATION MATRIX TAKING VECTOR IN LAB FRAME TO SAMPLE FRAME
% NECESSARY TO GET THE ORIENTATION OF CRYSTALS WITH RESPECT TO SAMPLE FRAME
RLab2Sam    = eye(3,3);

% SMOOTHING METHOD
% 1: Discrete Delta
% 2: Smoothed DiscreteDelta (slow)
% 3: AggregateFunction (medium)
SmoothingMethod = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load workspace for fundamental region.
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname)

% Load MIDAS results
pfname  = fullfile(pname, fname);
Grains  = parseGrainData(pfname, ws.frmesh.symmetries);
numpts  = length(Grains);
wts     = ones(1, numpts);

% CONVERT [R] TO SAMPLE COORDINATE SYSTEM
for i = 1:1:numpts
    RMats(:,:,i)   =  RLab2Sam*Grains(i).RMat;
end
quat    = QuatOfRMat(RMats);
rod     = ToFundamentalRegion(quat, ws.frmesh.symmetries);

% ODF GENERATION FROM DATA
if SmoothingMethod == 1
    tic
    [elem, ecrd] = MeshCoordinates(ws.frmesh, rod);
    odf = DiscreteDelta(ws.frmesh, ws.frmesh.l2ip, elem, ecrd, wts);
    odf = odf./MeanValue(odf, ws.frmesh.l2ip);
    t = toc;
    disp(['Time for DiscreteDelta:  ', num2str(t)]);
    
    PlotFR(ws.frmesh, odf, 'ShowMesh', 'on')
elseif SmoothingMethod == 2
    % Smoothed DiscreteDelta (slow)
    tic
    gqrule = QRuleGlobal(ws.frmesh, ws.frmesh.qrule, @RodMetric);
    
    [elem, ecrd] = MeshCoordinates(ws.frmesh, rod);
    odf = DiscreteDelta(ws.frmesh, ws.frmesh.l2ip, elem, ecrd, wts);
    odf = odf./MeanValue(odf, ws.frmesh.l2ip);
    
    odf = SmoothFunction(odf, ws.frmesh, ...
        ws.frmesh.qrule.pts, gqrule, ...
        @RodGaussian, std_sm, ws.frmesh.symmetries);
    odf = odf ./MeanValue(odf, ws.frmesh.l2ip);
    t = toc;
    disp(['Time for DiscreteDelta:  ', num2str(t)]);
    
    PlotFR(ws.frmesh, odf, 'ShowMesh', 'on')
elseif SmoothingMethod == 3
    % AggregateFunction (medium)
    pts   = ws.frmesh.crd(:, 1:ws.frmesh.numind);
    agg   = rod;
    wts   = ones(1, numpts);
    PointFun = @RodGaussian;
    
    tic
    odf = AggregateFunction(pts, agg, wts, PointFun, std_agg, ws.frmesh.symmetries);
    odf = odf./MeanValue(odf, ws.frmesh.l2ip);
    t = toc;
    disp(['Time for AggregateFunction:  ', num2str(t)]);
    
    PlotFR(ws.frmesh, odf, 'ShowMesh', 'on')
else
    disp('no smooothing method specified')
    PlotFRPerimeter('cubic')
    scatter3(rod(1,:), rod(2,:), rod(3,:), 30, 'b', 'filled')
    axis equal tight off
end