function hexp = HexPolytope()
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
piby12 = pi/12;
ztop = tan(piby12);
vnrm = sec(piby12);
%
%  Do top, then multiply by -1 to get the bottom.
%
vangs = piby12 + (0:11)*pi/6;
vtop  = [vnrm*cos(vangs); vnrm*sin(vangs); repmat(ztop, [1 12])];
vbot  = vtop; vbot(3, :) = -1*vbot(3, :);
%
vertices  = [vtop vbot];
%
tface  = 1:12;      tnormal = [0 0  1]; trhs = ztop; % top
bface  = 24:-1:13;  bnormal = [0 0 -1]; brhs = ztop; % bottom
%
lfaces  = [1:12; [12 1:11]; [24 13:23]; 13:24]';      % lateral faces
angs    = (0:11)*pi/6;
lnormal = [cos(angs); sin(angs); zeros(1, 12)]';
lrhs    = ones(12, 1);
%
hfaces = {tface, bface};
for i=1:12
  hfaces = {hfaces{:}, lfaces(i, :)};
end
%
matrix = [tnormal; bnormal; lnormal];
rhs    = [trhs; brhs; lrhs];
%
hexp = PolytopeStructure(matrix, rhs, vertices, hfaces);
