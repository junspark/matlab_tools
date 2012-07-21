function fsm = SmoothFunction(f, mesh, qpts, gqrule, SmoothFun, varargin)
% SmoothFunction - Smooth by convolution.
%   
%   USAGE:
%
%   fsm = SmoothFunction(f, mesh, gqrule, @SmoothFun)
%   fsm = SmoothFunction(f, mesh, gqrule, @SmoothFun, arg1, ...)
%
%   INPUT:
%
%   f      is an n-vector, 
%          the values of a function on the mesh at the
%          set of independent nodal points
%   mesh   is a MeshStructure,
%          for any region
%   qpts   is d x n,
%          the barycentric coordinates of the quadrature points
%   gqrule is a QRuleStructure,
%          the global quadrature rule for the mesh
%
%   SmoothFun is a function handle,
%             which defines the function to use for
%             smoothing; it has the form:
%
%             SmoothFun(center, points, ...)
%
%             center is m x 1, 
%                    the center of the distribution
%             points is m x n, 
%                    the list of points to evaluate
%
%   Other arguments are passed to `SmoothFun'
%
%   OUTPUT:
%
%   fsm is 1 x n, 
%      the nodal point values of the smoothed function at
%      the set of independent nodal points
%
%   NOTES:
% 
%   *  The function is smoothed by convolving a fixed distribution
%      (e.g. Gaussian) with the given function.  The result at 
%      each nodal point value is the integral of the input function 
%      times the smoothing function centered at that nodal point.
%
%   *  The result is not normalized.  
%      
%
f = f(:)';
CheckMesh(mesh, f);
%
crd = mesh.crd; dim = size(crd, 1);
con = mesh.con; nel = size(con, 2); nnpe = size(con, 1);
nqp = size(qpts, 2);
%
numind = length(f);
%
gqps = gqrule.pts; % Global quadrature points
gqpw = gqrule.wts; % and weights
%
els   = repmat(1:nel, [nqp 1]);
ecrds = repmat(qpts, [1 nel]);
fgqp  = EvalMeshFunc(mesh, f, els(:)', ecrds);
%
fwts = fgqp .* gqpw;
%
fsm = zeros(1, numind);
for i=1:numind
  disp(['node: ', num2str(i)])
  sfqp   = feval(SmoothFun, crd(:, i), gqps, varargin{:});
  fsm(i) = dot(fwts, sfqp);
end
%
