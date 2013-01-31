function etype = ElementTypeStruct(typename)
% ElementTypeStruct - Element type structure
%
%   USAGE:
%
%   etype = ElementTypeStruct(typename)
%   
%   INPUT:
%
%   typename is a string
%            Choices are:
%            (0-D)   'point'
%            (1-D)   'lines:2', 'lines:3'
%            (2-D)   'triangles:3', 'triangles:6', 'quads:4'
%            (3-D)   'tets:4', 'tets:10', 'bricks:8'
%
%   OUTPUT:
%
%   NOTES:
%
%   * Standard ordering of all elements is dictionary(z, y, x)  [z is sorted first]
%     i.e. x varies fastest, then y, then z 
%
%   * Need to do 'quads:8', 'quads:9', 'bricks:20'
%
if nargin < 1
  error('typename not specified')
end
%
tname = lower(typename);
switch tname
  %
  %  0-D
  %
 case 'point'
  dimension      = 0;
  is_rectangular = 1;
  is_simplicial  = 1;
  surftype = '';
  elsurfs  = [];
  %
  %  1-D
  %
 case 'lines:2'
  dimension      = 1;
  is_rectangular = 1;
  is_simplicial  = 1;
  surftype = 'point';
  elsurfs  = [1 2];
 case 'lines:3'
  dimension      = 1;
  is_rectangular = 1;
  is_simplicial  = 1;
  surftype = 'point';
  elsurfs  = [1 3];
  %
  %  2-D
  %
 case 'triangles:3'
  dimension      = 2;
  is_rectangular = 0;
  is_simplicial  = 1;
  surftype = 'lines:2';
  elsurfs  = [1 2; 3 1; 2 3]';
 case 'triangles:6'
  dimension      = 2;
  is_rectangular = 0;
  is_simplicial  = 1;
  surftype = 'lines:3';
  elsurfs  = [1 2 3; 6 4 1; 3 5 6]';
 case 'quads:4'
  dimension      = 2;
  is_rectangular = 1;
  is_simplicial  = 0;
  surftype = 'lines:2';
  elsurfs  = [1 2; 3 1; 2 4; 4 3]';
  %
  %  3-D
  %  surfaces ordered to produce inner normal
  %
 case 'tets:4'
  dimension      = 3;
  is_rectangular = 0;
  is_simplicial  = 1;
  surftype = 'triangles:3';
  elsurfs  = [1 2 3; 1 3 4; 1 4 2; 2 4 3]';
 case 'tets:10'  % femevps
  dimension      = 3;
  is_rectangular = 0;
  is_simplicial  = 1;
  surftype = 'triangles:6';
  elsurfs  = [1 2 3 4 5 6; 1 4 6 7 9 10; 1 7 10 2 8 3; 3 8 10 5 9 6]';
 case 'bricks:8'  
  dimension      = 3;
  is_rectangular = 1;
  is_simplicial  = 0;
  surftype = 'quads:4';
  elsurfs  = [1 2 3 4; 2 1 6 5; 1 3 5 7; 4 2 8 6; 3 4 7 8; 7 8 5 6]'; 
  %
  %  case 'bricks:20'  TO BE DONE
  %
 otherwise
  error(sprintf('no such element type:  %s', typename))
end
%
etype.name     = tname;
etype.surfs    = elsurfs;
etype.surftype = surftype;
%
etype.dimension      = dimension;
etype.is_rectangular = is_rectangular;
etype.is_simplicial  = is_simplicial;
