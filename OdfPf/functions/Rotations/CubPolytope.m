function cubp = CubPolytope()
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
     x  x  z;  x  z  x;  z  x  x;
    -x  x  z; -x  z  x; -z  x  x;
    -x -x  z; -x -z  x; -z -x  x;
     x -x  z;  x -z  x;  z -x  x;
     x  x -z;  x  z -x;  z  x -x;
    -x  x -z; -x  z -x; -z  x -x;
    -x -x -z; -x -z -x; -z -x -x;
     x -x -z;  x -z -x;  z -x -x;
	   ]';
%
%  Form polygons for faces.
%
f100 = {...
    [ 1     2    11    10    22    23    14    13], ...
    [ 4     5     8     7    19    20    17    16], ...
    [ 1     3     6     4    16    18    15    13], ...
    [ 7     9    12    10    22    24    21    19], ...
    [ 2     3     6     5     8     9    12    11], ...
    [14    15    18    17    20    21    24    23] ...
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
