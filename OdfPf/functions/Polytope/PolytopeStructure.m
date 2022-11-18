function pstruct = PolytopeStructure(matrix, rhs, vertices, faces)
% PolytopeStructure - Structure for polytope.
%
%   USAGE:
%
%   pstruct = PolytopeStructure(matrix, rhs)
%   pstruct = PolytopeStructure(matrix, rhs, vertices)
%   pstruct = PolytopeStructure(matrix, rhs, vertices, faces)
%
%   INPUT:
%
%   matrix   is a m x d real array:
%            the matrix of the linear inequalities defining
%            the polytope (A, where Ax <= b)
%   rhs      is a m x 1 real array:
%            the right hand side of the linear inequalities defining
%            the polytope (b, where Ax <= b)
%   vertices is a d x n real array:  (optional)
%            the list of n vertices in R^d
%   faces    is a 1 x f cell array:  (optional)
%            each cell contains a list of vertex numbers
%            for a given face (3D polytopes only)
%
%   OUTPUT:
%
%   pstruct is a Polytopestructure,
%           it consists of two primary fields, .matrix and
%           .rhs; optionally, it can have two more fields,
%           .vertices and .faces; the optional fields
%           are not computed here, and they are primarily
%           used by the graphics functions
%   
pstruct.matrix   = matrix;
pstruct.rhs      = rhs;
%
if (nargin >= 3)
  pstruct.vertices = vertices;
else
  pstruct.vertices = [];
end % if (nargin >= 3)
%
if (nargin >= 4)
  pstruct.faces = faces;
else
  pstruct.faces = [];
end % if (nargin >= 4)
