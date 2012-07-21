function newmesh = ReduceMesh(mesh, sym, tol)
% ReduceMesh - Find equivalence array for a Rodrigues mesh.
%   
%   USAGE:
%
%   newmesh = ReduceMesh(mesh, sym)
%   newmesh = ReduceMesh(mesh, sym, tol)
%
%   INPUT:
%
%   mesh is a MeshStructure on the fundamental region
%        with (possibly) an empty equivalence array
%   sym  is 4 x m, 
%        the symmetry group in quaternions
%   tol  is 1 x 1, (optional, default: 10^-14) 
%        the tolerance used for comparing equivalent rotations; 
%
%   OUTPUT:
%
%   newmesh is a MeshStructure on the fundamental region;
%           it has the same  and connectivity, but
%           reordered; the equivalence array is created
%           and the nodes are numbered so that all 
%           the independent nodes precede the dependent nodes.
%
crd = mesh.crd;
con = mesh.con;
%
n = size(crd, 2);
%
if (n <= 1) % kind of silly ...
  newmesh = mesh;
  newmesh.eqv = [];
  return
end
% 
if (nargin < 2)
  error('require symmetry group (arg 2)');
end
%
if (nargin < 3)
  tol = 1.0e-14;  % default value of tolerance
end
%
%--------------------Real work starts here------------------------------
%
%  Find misorientation with four other orientations,
%  and sort uniquely to determine equivalences.
%
angs = [0 exp(1./(2:4))];
%
R3 =[0.6685   -0.3853    0.6362;
     0.6649    0.6928   -0.2791;
     -0.3332    0.6095    0.7193;
     ];
a1 = [1 0.5 0.2]';
a2 = R3*a1;
a3 = R3*a2;
%
e123 = [a1 a2 a3];
a111 = [1 1 1]';
axxx = [a111 e123];
%
fourpts = QuatOfAngleAxis(angs, axxx);
% 
qcrd = QuatOfRod(crd);
%
miso = zeros(4, n);
for i=1:4
  miso(i, :) = Misorientation(fourpts(:, i), qcrd, sym);
end
%
[u, ord, iord] = UniqueVectors(miso, tol);
nind  = length(ord);
%
%  Mark independent and dependent nodes.
%
eqvnp      = ones(1, n);
eqvnp(ord) = 0;
indnp      = (eqvnp == 0); % independent nodal points
depnp      = ~indnp;
%
%  Generate equivalence and numbering for dependent nodes.
%
eqvnum     = cumsum(eqvnp);
eqv(1, :)  = eqvnum(depnp) + nind;
eqv(2, :)  = iord(depnp);
%
%  Complete the numbering.
%
number(indnp) = iord(indnp);
number(depnp) = eqv(1, :);
%
[tmp, inumber] = sort(number);
newcrd  = crd(:, inumber);
newcon  = number(con);
newmesh = MeshStructure(newcrd, newcon, eqv);

