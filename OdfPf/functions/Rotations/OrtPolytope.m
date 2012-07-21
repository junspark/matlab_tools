function ortp = OrtPolytope()
% OrtPolytope - Polytope for orthorhombic symmetry.
%   
%   USAGE:
%
%   ortp = OrtPolytope
%
%   INPUT:  none
%
%   OUTPUT:
%
%   ortp is a PolytopeStructure:
%        it gives the polytope for the orthorhombic 
%        fundamental region including the vertex
%        list and the faces component (for plotting)
%

%  
%  Compute faces (constraints).
%
matrix = [eye(3); -eye(3)];
rhs    = ones(6, 1);
%
%  Compute Vertices.
%
vertices = [-1 -1 -1;  1 -1 -1; 1  1 -1; -1 1 -1;...
	    -1 -1  1;  1 -1  1; 1  1  1; -1 1  1]';
%
faces = {...
    [4 3 2 1], [5 6 7 8], [2 3 7 6], ...
    [5 8 4 1], [1 2 6 5], [7 8 4 3]  ...
	};
%
ortp = PolytopeStructure(matrix, rhs, vertices, faces);
