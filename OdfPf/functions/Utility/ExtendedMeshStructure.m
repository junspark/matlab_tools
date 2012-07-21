function extMesh = ExtendedMeshStructure(meshStructure, quadRule, varargin)
%
%  USAGE
%    extMesh = ExtendedMeshStructure(meshStructure, quadRule, rigthSymGroup, {leftSymGroup})
%
%       crd: the mesh crd's
%       con: the connectivity
%       eqv: the equivalence array (x-stal symmetry only)
%     h1sip: the H1 semi inner-product matrix
%      l2ip: the L2 innerproduct matrix
%     qrule: the reference quadrature rule
%      rsym: the right(crystal) symmetry group in quaternion representation
%      lsym: the left (sample/additional lattice) symmetry group in
%            quaternion representation 
% NOTES
%     *) Specialized to simplex elements for now...

% ------- defaults -------
lsymGroupFlag = 0;
lsymGroup = [];
rsymGroupFlag = 0;
rsymGroup = [];
forceSphFlag = 0;
% ------- defaults -------

% Process options
if nargin > 2
  for i_arg = 1:2:length(varargin)
    if strcmp(varargin{i_arg}, 'LeftSym')
      lsymGroupFlag = 1;
      lsymGroup = varargin{i_arg + 1};
    elseif  strcmp(varargin{i_arg}, 'RightSym')
      rsymGroupFlag  = 1;
      rsymGroup = varargin{i_arg + 1};
    elseif strcmp(varargin{i_arg}, 'UseSph')
      if strcmp(varargin{i_arg + 1}, 'on')
        forceSphFlag = 1;
      end
    end
  end
end
dim  = size(meshStructure.crd, 1);
edim = size(meshStructure.con, 1);
%
fprintf('\nChecking for inside-out elements...');
meshStructure = FixConnectivity(meshStructure);
%
if edim == dim + 1 & ~forceSphFlag
  %
  fprintf('\nForming L2IP Matrix on the Rodrigues space...');
  tic
  l2ip = L2IPMatrix(meshStructure, quadRule, @RodMetric);
  fprintf('done\n');
  toc
  %
  fprintf('\nForming H1SIP Matrix on the Rodrigues space...');
  tic
  [grad, metricGij] = NpQpGradMatrix(meshStructure, quadRule.pts);
  h1sip = H1SIPMatrix(meshStructure, quadRule, grad, metricGij);
  fprintf('done\n');
  toc
  %
elseif edim == dim + 1 & forceSphFlag
  qmeshStructure = meshStructure;
  qmeshStructure.crd = QuatOfRod(meshStructure.crd);
  %
  fprintf('\nForming L2IP Matrix on the sphere...');
  tic
  l2ip = SphL2IP(qmeshStructure, quadRule);
  fprintf('done\n');
  toc
  %
  fprintf('\nForming H1SIP Matrix on the sphere...');
  tic
  [grad, metricGij] = NpQpGradMatrixSph(qmeshStructure, quadRule.pts);
  h1sip = H1SIPMatrixSph(qmeshStructure, quadRule, grad, metricGij);
  fprintf('done\n');
  toc
  %
elseif edim == dim
  %
  fprintf('\nForming L2IP Matrix on the sphere...');
  tic
  l2ip = SphL2IP(meshStructure, quadRule);
  fprintf('done\n');
  toc
  %
  fprintf('\nForming H1SIP Matrix on the sphere...');
  tic
  [grad, metricGij] = NpQpGradMatrixSph(meshStructure, quadRule.pts);
  h1sip = H1SIPMatrixSph(meshStructure, quadRule, grad, metricGij);
  fprintf('done\n');
  toc
  %
end

extMesh       = meshStructure;
extMesh.l2ip  = l2ip;
extMesh.h1sip = h1sip;
extMesh.qrule = quadRule;
extMesh.rsym  = rsymGroup;
extMesh.lsym  = lsymGroup;
%
