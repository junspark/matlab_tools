function mesh = SphBaseMesh(dim, varargin)
% SphBaseMesh - Generate base mesh for spheres.
%
%   USAGE:
%
%   mesh = SphBaseMesh(dim)
%   mesh = SphBaseMesh(dim, 'param', 'value')
%
%   INPUT:
%
%   dim is a positive integer,
%       the dimension of the sphere (2 for the usual sphere S^2)
%
%   These arguments can be followed by a list of
%   parameter/value pairs which specify keyword options.
%   Available options include:
%
%   'Hemisphere'  'on'|'off'
%                 to mesh only the upper hemisphere
%   
%
%   OUTPUT:
%
%   mesh is a MeshStructure,
%        on the sphere of the specified dimension
%
%   NOTES:
%
%   * For S^2, the normals may be mixed, some outer, some inner.  
%     This needs to be fixed.
%

%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'Hemisphere', 'off' ...
       };
%
opts = OptArgs(optcell, varargin);
%
HEMI = OnOrOff(opts.Hemisphere);
%
%
%-------------------- Form Base Mesh
%
n = dim;
%
%  Set up points.  If hemisphere, ignore last one.
%
caxes  = diag(1:(n+1));
pts    = [zeros(n+1, 1), caxes, -caxes];
conlen = 2^(n+1);
if (HEMI)
  pts    = pts(:, 1:(end - 1));
  conlen = 2^n;
end
%
%  Form connectivity by using delaunayn on a fully 
%  dimensional simplex, and then remove the origin.
%
con  = delaunayn(pts')' - 1;
con  = reshape(con(con > 0), [n+1, conlen]);
mesh = MeshStructure(UnitVector(pts(:, 2:end)), con);
