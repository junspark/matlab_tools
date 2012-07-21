function odf = ProcessAggregate(ws, quat, wts, method, varargin)
%  ProcessAggregate - Build an ODF from an aggregate
%
%   VERSION:  $Id: ProcessAggregate.m 169 2010-02-17 19:23:30Z boyce $
%   STATUS:  in development
%
%   USAGE:
%
%   odf = ProcessAggregate(ws, quat, wts, method, <options>)
%
%   INPUT:
%
%   ws is a workspace 
%         with fundamental region mesh and symmetries
%   quat is 4 x n
%         array of quaternions (the aggregate)
%   wts  is 1 x n (or empty)
%         weights to use; empty weights means uniform
%   method is a string
%         either 'DiscreteDelta' or 'AggregateFunction' (capitalization does not matter)
%
%   <options> is a sequence of parameter, value pairs
%             Available options are listed below.  Default values,
%             if set, are shown in brackets.
%
%       'GaussianRadius'    scalar [5 * pi/180]
%                           smoothing angle for aggregate function
%       'Verbose'           {true} | false
%
%   OUTPUT:
%
%   odf is a k-vector  (where k is number of independent nodes)
%          the resulting ODF 
%
%   NOTES:
%
%   *  This routine simplifies the interface to DiscreteDelta or AggregateFunction
%      for the most common cases
%
MyName = mfilename;
%
%-------------------- Defaults and Keyword Arguments
%
ToRadians = pi/180;
optcell = {...
    'GaussianRadius',  ToRadians*5, ...
    'Verbose',         true ...
       };
%
options = OptArgs(optcell, varargin);

%
rod = ToFundamentalRegion(quat, ws.frmesh.symmetries);
if isempty(wts)
  wts = ones(1, size(quat, 2));
end

switch lower(method)
  %
 case 'discretedelta'
  %
  if options.Verbose
    disp 'using discrete delta'
  end
  %
  [elem, ecrd] = MeshCoordinates(ws.frmesh, rod);
  odf = DiscreteDelta(ws.frmesh, ws.frmesh.l2ip, elem, ecrd, wts);
  %
 case 'aggregatefunction'
  %
  if options.Verbose
    disp 'using aggregate function'
  end
  %
  odf = AggregateFunction(ws.frmesh.crd(:, 1:ws.frmesh.numind), rod, wts, ...
			  @RodGaussian, options.GaussianRadius, ws.frmesh.symmetries);
  %
 case 'positivediscretedelta'
  %
  if options.Verbose
    disp 'using positive discrete delta'
  end
  ncrd = size(ws.frmesh.crd, 2);
  %
  npts   = length(elem);
  [m, n] = size(ecrd);
  %
  mn = m*n;
  i = ws.frmesh.con(:, elem); i = i(:);
  j = ones(mn, 1);
  w = repmat(wts(:)', [m 1]); w = w(:);
  v = ecrd(:).*w;
  %
  rhs = sparse(i, j, v, ncrd, 1);
  rhs = full(rhs);
  %  
  %  Now reduce according to equivalence.
  %
  rhs = EqvReduce(rhs', ws.frmesh.eqv);
  %
  objfun = @(x) pddObjective(ws.frmesh.l2ip, rhs(:), x);
  x0 = DiscreteDelta(ws.frmesh, ws.frmesh.l2ip, elem, ecrd, wts);
  x0(x0 < 0) = 0;
  lb = zeros(length(rhs), 1);
  opts = optimset('Display', 'iter', 'GradObj', 'on');
  odf = fmincon(objfun, x0, [], [], [], [], lb, [], [], opts);
  %
 otherwise
  error(sprintf('no such method:  %s', method))
end
%
%  Normalize to mean value one (multiples of uniform)
%  * make a column vector on return
%
odf = odf(:) ./ MeanValue(odf, ws.frmesh.l2ip);


function [f, g] = pddObjective(mat, rhs, x)
% PDDOBJECTIVE - 
%   
mx = mat*x - rhs;
f = mx'*mx;
g = 2*mat*mx;
