function hexp = TrigPolytope()
% HexPolytope - Polytope for hexagonal fundamental region.
%   
%   USAGE:
%
%   hexp = HexPolytope
%
%   INPUT:  none
%
%   OUTPUT:
%
%   hexp is a PolytopeStructure:
%        the polytope of the hexgonal fundamental region;
%        it includes the vertex list and the list of
%        polygonal faces
%
piby6 = pi/6;
ztop = tan(piby6);
vnrm = sec(piby6);
%
%  Do top, then multiply by -1 to get the bottom.
%
vangs = piby6 + (0:5)*pi/3;
vtop  = [vnrm*cos(vangs); vnrm*sin(vangs); repmat(ztop, [1 6])];
vbot  = vtop; vbot(3, :) = -1*vbot(3, :);
%
vertices  = [vtop vbot];
%
tface  = 1:6;      tnormal = [0 0  1]; trhs = ztop; % top
bface  = 12:-1:7;  bnormal = [0 0 -1]; brhs = ztop; % bottom
%
lfaces  = [1:6; [6 1:5]; [12 7:11]; 7:12]';      % lateral faces
angs    = (0:5)*pi/3;
lnormal = [cos(angs); sin(angs); zeros(1, 6)]';
lrhs    = ones(6, 1);
%
hfaces = {tface, bface};
for i=1:6
  hfaces = {hfaces{:}, lfaces(i, :)};
end
%
matrix = [tnormal; bnormal; lnormal];
rhs    = [trhs; brhs; lrhs];
%
hexp = PolytopeStructure(matrix, rhs, vertices, hfaces);
