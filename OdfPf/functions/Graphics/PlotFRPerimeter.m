function perim =  PlotFRPerimeter(symtype)
% PlotFRPerimeter - Plot perimeter of fundamental region.
%   
%   USAGE:
%
%   perim =  PlotFRPerimeter(symtype)
%
%   INPUT:
%
%   symtype is a string,
%           indicating the fundamental region type; possible
%           values are:  'cubic'|'hexagonal'|'orthorhombic'; 
%           if the symmetry type is not recognized, a warning 
%           is issued
%
%   OUTPUT:
%
%   perim is a column vector of line handles:
%            it contains handles to the lines outlining each face,
%            as returned by the matlab `plot3' function
%
%   NOTES:
%
%   * sets `hold' to on
%
if (strcmpi(symtype, 'cubic'))
  ptope = CubPolytope;
elseif (strcmpi(symtype, 'hexagonal'))
  ptope = HexPolytope;
elseif (strcmpi(symtype, 'orthorhombic'))
  ptope = OrtPolytope;
else
  wid = 'PlotFRPerimeter:symmetry';
  msg = 'symmetry type not recognized';
  warning(wid, msg)
  return
end
%
verts  = ptope.vertices;
faces  = ptope.faces;
nfaces = length(faces);
%
perim = [];
for i=1:nfaces
  list = faces{i};
  pts = verts(:, [list list(1)]);
  prm = plot3(pts(1, :), pts(2, :), pts(3, :), 'ko-');
  perim = [perim; prm];
  hold on
end
