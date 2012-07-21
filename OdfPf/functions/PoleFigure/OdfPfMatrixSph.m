function opm = OdfPfMatrixSph(hkl, mesh, sym, pts, div, varargin)
% ODFPFMATRIX - Find matrix relating ODF to pole figure.
%
%   opm = OdfPfMatrix(hkl, mesh, sym, pts, div, invpf)
%   opm = OdfPfMatrix(hkl, mesh, sym, pts, div)
%
%   hkl   is 3 x 1, the pole direction
%   mesh  is a MeshStructure on orientation space
%   sym   is 4 x s, the symmetry group in quaternions
%   pts   is 3 x p, a list of p points on the sphere
%   div   is 1 x 1, the number of divisions per fiber
%   invpf is 1 x 1, (optional) a nonzero value flags causes
%                   computation of inverse pole figure matrix
%
%   opm is p x n, (sparse) the matrix which takes nodal point values
%                 on the fundamental region to pole figure/inverse
%                 pole figure values on the sphere
%
%   Notes:
%
%   *)  opm acts on the "reduced" set of nodes, the set of nodes
%       in which equivalent nodes are combined into a single
%       degree of freedom
%   *)  this routine is very memory intensive; you may need to
%       build the matrix in pieces; if you have a lot of points
%       and a large value of `div', then you should break the points
%       into smaller groups
%
%
crd = mesh.crd;
con = mesh.con;
eqv = mesh.eqv;

num_tnp = size(crd, 2);      % total number of nodes
num_pts = size(pts, 2);      % number of pole figure points

opm = sparse(div*num_pts, num_tnp);

% Defaults
invfib = 0;
userLeftSymFlag = 0;

% Process cases
if nargin > 5
  for i_arg = 1:2:length(varargin)
    if strcmp(varargin{i_arg}, 'InversePF')
      if strcmp(varargin{i_arg + 1}, 'on')
        invfib = 1;
      end
    elseif strcmp(varargin{i_arg}, 'LeftSym')
      if strcmp(varargin{i_arg + 1}, 'axisymmetry')
        userLeftSymFlag = 1;
        leftSymTag = varargin{i_arg + 1};
      end
    elseif strcmp(varargin{i_arg}, 'symmAxis')
      symmAxis = varargin{i_arg + 1};
    end
  end
end

% Find the fiber of {hkl}||pts
fib = FiberOfPoint(pts, hkl, div, varargin{:}, 'QuatOut', 'on');

% Check for the special case of axisymmetry
if userLeftSymFlag
  fib = ToAxisymmetricFR(reshape(fib, [4, div*num_pts]), sym, symmAxis);
end

% Get the odfpf matrix from the sph interpolation matrix
opm = InterpMatrixOfSphMesh(mesh, reshape(fib, [size(fib, 1), div*num_pts]));
opm = reshape(opm, [div, num_pts*num_tnp]);
opm = sum(opm, 1)/div;
opm = reshape(opm, [num_pts, num_tnp]);

% Reduce by equivalence...
opm = EqvReduce(opm, eqv);

