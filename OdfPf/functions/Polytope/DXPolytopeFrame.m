function h = DXPolytopeFrame(ptope, prefix)
% DXPOLYTOPEFRAME - Write polytope frame to DX file.
%   
%   h = DXPolytopeFrame(p, prefix)
%
h = {};
%
%  positions
%
v  = ptope.vertices;
vl = size(v, 2);
fn = [prefix, '-vertices.dat'];
h = {h{:}, ...
     'object "vertex list" class array type float', ...
     ['rank 1 shape 3 items ', num2str(vl)], ...
     ['ascii data file ', fn] ...
     };
WriteDataFile(fn, {}, v');
%
%  connections
%
f = ptope.faces;  % cell array of vertex lists
con = [];
for i=1:length(f)
  facei = f{i};
  li = length(facei);
  con = [con, [facei(1:li); [facei(2:li), facei(1)]]];
end
%
fn = [prefix, '-edges.dat'];
lc = size(con, 2);
h = {h{:}, ...
     'object "edge list"   class array    type int', ...
     ['rank 1 shape 2 items ', num2str(lc)], ...
     ['ascii data file ', fn], ...
     'attribute "element type" string "lines"', ...
     'attribute "ref" string "positions"'
     };
%
WriteDataFile(fn, {}, (con - 1)', 'integer');
%
%  Make DX field.
%
h = {h{:}, ...
     ['object "', [prefix, '-polytope'],'" class field'],...
     '   component "positions" "vertex list"',  ...
     '   component "connections" "edge list"', ...
     'end' ...
     };
%
WriteCellStrings([prefix, '.dx'], h);
%          0          1
%          1         10
%         10          9
%          9         21
%         21         22
%         22         13
%         13         12
%         12          0
%          3          4
%          4          7
%          7          6
%          6         18
%         18         19
%         19         16
%         16         15
%         15          3
%          0          2
%          2          5
%          5          3
%          3         15
%         15         17
%         17         14
%         14         12
%         12          0
%          6          8
%          8         11
%         11          9
%          9         21
%         21         23
%         23         20
%         20         18
%         18          6
%          1          2
%          2          5
%          5          4
%          4          7
%          7          8
%          8         11
%         11         10
%         10          1
%         13         14
%         14         17
%         17         16
%         16         19
%         19         20
%         20         23
%         23         22
%         22         13
%          0          1
%          1          2
%          2          0
%          3          4
%          4          5
%          5          3
%          6          7
%          7          8
%          8          6
%          9         10
%         10         11
%         11          9
%         12         13
%         13         14
%         14         12
%         15         16
%         16         17
%         17         15
%         18         19
%         19         20
%         20         18
%         21         22
%         22         23
%         23         21
%    4.142135623730950e-01    4.142135623730950e-01    1.715728752538099e-01
%    4.142135623730950e-01    1.715728752538099e-01    4.142135623730950e-01
%    1.715728752538099e-01    4.142135623730950e-01    4.142135623730950e-01
%   -4.142135623730950e-01    4.142135623730950e-01    1.715728752538099e-01
%   -4.142135623730950e-01    1.715728752538099e-01    4.142135623730950e-01
%   -1.715728752538099e-01    4.142135623730950e-01    4.142135623730950e-01
%   -4.142135623730950e-01   -4.142135623730950e-01    1.715728752538099e-01
%   -4.142135623730950e-01   -1.715728752538099e-01    4.142135623730950e-01
%   -1.715728752538099e-01   -4.142135623730950e-01    4.142135623730950e-01
%    4.142135623730950e-01   -4.142135623730950e-01    1.715728752538099e-01
%    4.142135623730950e-01   -1.715728752538099e-01    4.142135623730950e-01
%    1.715728752538099e-01   -4.142135623730950e-01    4.142135623730950e-01
%    4.142135623730950e-01    4.142135623730950e-01   -1.715728752538099e-01
%    4.142135623730950e-01    1.715728752538099e-01   -4.142135623730950e-01
%    1.715728752538099e-01    4.142135623730950e-01   -4.142135623730950e-01
%   -4.142135623730950e-01    4.142135623730950e-01   -1.715728752538099e-01
%   -4.142135623730950e-01    1.715728752538099e-01   -4.142135623730950e-01
%   -1.715728752538099e-01    4.142135623730950e-01   -4.142135623730950e-01
%   -4.142135623730950e-01   -4.142135623730950e-01   -1.715728752538099e-01
%   -4.142135623730950e-01   -1.715728752538099e-01   -4.142135623730950e-01
%   -1.715728752538099e-01   -4.142135623730950e-01   -4.142135623730950e-01
%    4.142135623730950e-01   -4.142135623730950e-01   -1.715728752538099e-01
%    4.142135623730950e-01   -1.715728752538099e-01   -4.142135623730950e-01
%    1.715728752538099e-01   -4.142135623730950e-01   -4.142135623730950e-01
% object "vertex list" class array type float
% rank 1 shape 3 items 24
% ascii data file cfr-vertices.dat
% object "edge list"   class array    type int
% rank 1 shape 2 items 72
% ascii data file cfr-edges.dat
% attribute "element type" string "lines"
% attribute "ref" string "positions"
% object "cfr-polytope" class field
%    component "positions" "vertex list"
%    component "connections" "edge list"
% end
