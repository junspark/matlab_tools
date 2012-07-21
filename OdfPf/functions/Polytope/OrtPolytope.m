function cubp = OrtPolytope()
% CubPolytope - Polytope for cubic fundamental region.
%   
%   USAGE:
%
%   cubp = CubPolytope
%
%   INPUT:  none
%
%   OUTPUT:
%
%   cubp is a PolytopeStructure:
%        it gives the polytope for the cubic 
%        fundamental region including the vertex
%        list and the faces component (for plotting)
%

%  
%  Compute faces (constraints).
%
b1 = tan(pi/8);
b2 = tan(pi/6);
%
n111   = UnitVector([1 1 1; 1 -1 1; -1 1 1; -1 -1 1]');
matrix = [eye(3); n111'];
matrix = [matrix; -matrix];
%
rhs = [b1 b1 b1 b2 b2 b2 b2];
rhs = [rhs'; rhs'];
%
%  Compute Vertices.
%
x = b1; z = b1^2;
vertices = [...
     1  1  1;  1  1 -1;  1 -1 -1;  1 -1  1;
    -1 -1 -1; -1 -1  1; -1  1  1; -1  1 -1;
	   ]';
%
%  Form polygons for faces.
%
faces = {...
    [ 1     2     3     4], ...
    [ 1     2     3     4], ...
    [ 1     2     3     4], ...
    [ 1     2     3     4], ...
    [ 1     2     3     4], ...
    [ 1     2     3     4], ...
       };
    %
f111 = {[1 2 3], [4 5 6], [7 8 9], [10 11 12], ...
	[13 14 15], [16 17 18], [19 20 21], [22 23 24]};
%
faces = {f100{:}, f111{:}};
%
%  Form output structure.
%
cubp = PolytopeStructure(matrix, rhs, vertices, faces);
