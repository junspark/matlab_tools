function odf = CalcOdfFromAggregate(frmesh, quat, wts, varargin)
% CalcOdfFromAggregate - calculates odf from aggregate data
%
%   USAGE:
%
%   odf = CalcOdfFromAggregate(frmesh, quat, wts)
%   odf = CalcOdfFromAggregate(frmesh, quat, wts, varargin)
%
%   INPUT:
%
%   frmesh
%       standard odfpf frmesh structure. higher refinement level generally
%       means longer execution time.
%
%   quat
%       list of quaternions representing the crystallographic orientations
%       of the constiuents in the aggregate (4 x n)
%
%   wts
%       weights for the list of quaternions. if equal weight for each point
%       is desired use ones(n,1). (n x 1)
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'SmoothingMethod'  smoothing method; options below
%                      'DiscreteDelta'
%                      'SmoothedDiscreteDelta'  slow
%                      'AggregateFunction'      medium : default
%   'std'              standard deviation used for smoothing. this number
%                      has no effect if the 'DiscreteDelta' method is 
%                      chosen.
%   'PlotOdf'          plots the computed odf over rodrigues space.
%                      'on'     plots result
%                      'off'    plot suppressed (default)
%
%   OUTPUT:  
%
%   odf
%       orientation ditribution function of the aggregate.
%
%   NOTE:
%           This function is orginally from ODFPF example
%           "OdfFromAggregate"
%

% default options
optcell = {...
    'SmoothingMethod', 'AggregateFunction', ...
    'std', 5, ...
    'PlotOdf', 'off', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

% smoothing angle (rad)
std_smoothing   = deg2rad(opts.std);

% smoothing method
SmoothingMethod = opts.SmoothingMethod;

% number of points
numpts  = size(quat,1);           

% symmetry operator on the input quaternions to find the equivalent
% rodrigues vector in the fundamental region
rod  = ToFundamentalRegion(quat, frmesh.symmetries);

% ODF GENERATION FROM DATA
if strcmpi(SmoothingMethod, 'DiscreteDelta')
    tic
    [elem, ecrd] = MeshCoordinates(frmesh, rod);
    odf = DiscreteDelta(frmesh, frmesh.l2ip, elem, ecrd, wts);
    odf = odf./MeanValue(odf, frmesh.l2ip);
    t = toc;
    disp(['Time for DiscreteDelta:  ', num2str(t)]);
    
    if strcmpi(opts.PlotOdf, 'on')
        PlotFR(frmesh, odf, 'ShowMesh', 'on')
    end
elseif strcmpi(SmoothingMethod, 'SmoothedDiscreteDelta')
    % Smoothed DiscreteDelta (slow)
    tic
    gqrule = QRuleGlobal(frmesh, frmesh.qrule, @RodMetric);
    
    [elem, ecrd] = MeshCoordinates(frmesh, rod);
    odf = DiscreteDelta(frmesh, frmesh.l2ip, elem, ecrd, wts);
    odf = odf./MeanValue(odf, frmesh.l2ip);
    
    odf = SmoothFunction(odf, frmesh, ...
        frmesh.qrule.pts, gqrule, ...
        @RodGaussian, std_smoothing, frmesh.symmetries);
    odf = odf ./MeanValue(odf, frmesh.l2ip);
    t = toc;
    disp(['Time for DiscreteDelta:  ', num2str(t)]);
    
    if strcmpi(opts.PlotOdf, 'on')
        PlotFR(frmesh, odf, 'ShowMesh', 'on')
    end
elseif strcmpi(SmoothingMethod, 'AggregateFunction')
    % AggregateFunction (medium)
    pts   = frmesh.crd(:, 1:frmesh.numind);
    agg   = rod;
    PointFun = @RodGaussian;
    
    tic
    odf = AggregateFunction(pts, agg, wts, PointFun, std_smoothing, frmesh.symmetries);
    odf = odf./MeanValue(odf, frmesh.l2ip);
    t = toc;
    disp(['Time for AggregateFunction:  ', num2str(t)]);
    
    if strcmpi(opts.PlotOdf, 'on')
        PlotFR(frmesh, odf, 'ShowMesh', 'on')
    end
else
    disp('no smooothing method specified')
    
    if strcmpi(opts.PlotOdf, 'on')
        PlotFRPerimeter('cubic')
        scatter3(rod(1,:), rod(2,:), rod(3,:), 30, 'b', 'filled')
        axis equal tight off
    end
end