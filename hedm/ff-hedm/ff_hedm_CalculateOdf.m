%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ODF From Aggregate - MIDAS RESULTS
%  THIS IS DEPRECATED - MATLAB FUNCTION VERSION AVAIALBE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DegToRad = (pi/180);     % matlab `deg2rad' requires toolbox
wsname = 'wscub4x';      % workspace name
std_agg = DegToRad * 5;  % std deviation for AggregateFunction
std_sm  = DegToRad * 10; % std deviation for smoothing

pname_root  = '/home/beams/S1IDUSER/mnt/orthros/internal_aug18_midas/ff/hedm_resolution';
pname0      = fullfile(pname_root, 'ss_sam_ff3_Layer8_Analysis_Time_2018_09_11_12_27_14');

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
grain_map   = parseGrainData_OneLayer(pname0, CubSymmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', nan, ...
    'ComputeSelfMisoTable', false, ...
    'ComputeSelfDistTable', false, ...
    'OutputReflectionTable', false, ...
    'NumFrames', 1440, ...
    'Technique', 'ff-midas');

numpts  = grain_map.nGrains;
wts     = ones(1, numpts);
quat    = [grain_map.grains(:).quat];
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