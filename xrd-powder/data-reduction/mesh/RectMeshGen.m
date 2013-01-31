function mesh = RectMeshGen(ndiv, varargin)
% RectMeshGen - Mesh a rectangular region.
%   
%   USAGE:
%
%   mesh = RectMeshGen(ndiv)
%   mesh = RectMeshGen(ndiv, <options>)
%
%   INPUT:
%
%   ndiv is an n-vector of (integer)
%        It gives the subdivisions in each coordinate direction.  
%        n (number of dimensions) can be 1, 2 or 3.
%
%   <options> is a sequence of parameter, value pairs
%             Available options are listed below.  Default values,
%             if set, are shown in brackets.
%
%   'Transformation'     n x m  matrix
%                        Linear transformation applied to nodal points.
%                        Defaults to identity.
%   'Translation'        n x 1 vector
%                        Translation of nodal points, applied after transformation.
%                        Defaults to zero.
%   'ElementType'        string
%                        Defaults to 'bricks:8'.  See 'ElementTypeStruct'.
%
%   OUTPUT:
%
%   mesh is a MeshStructure 
%        It is for the transformed unit n-cube and uses the simplest rectangular
%        elements.
%
%   NOTES:
%
MyName = mfilename;
%
nd = length(ndiv);
if ( (nd ~= 1) &&  (nd ~= 2) &&  (nd ~= 3) )
  error('dimension must be 1, 2 or 3');
end
%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'Transformation',  eye(nd),      ...
    'Translation',     zeros(nd, 1), ...
    'ElementType',     'bricks:8'    ...
       };
%
options = OptArgs(optcell, varargin);
%
%  Construct meshes by increasing dimension.  This could 
%  be done more generally with a loop, but this is good 
%  enough for now.
%
%  Dimension 1.
%  -----------
%
n1   = ndiv(1);
dx   = 1.0/n1;
xcrd = 0:dx:1;
%
c1  = 1:n1;
c2  = c1 + 1;
con = [c1; c2];
%
if (nd == 1)
  crd  = xcrd;
  mesh = MeshStructure(crd, con, [], 'ElementType', 'lines:2');
  try
    mesh = ApplyTransformation(mesh, options.Transformation, ...
			       options.Translation);
  catch
    error(lasterr);
  end
  %
  mesh.ndiv = ndiv;
  return
end
%
%  Dimension 2.
%  -----------
%
n2   = ndiv(2);
dy   = 1.0/n2;
ycrd = 0:dy:1;
%
[xtmp, ytmp] = ndgrid(xcrd, ycrd);
xtmp = xtmp(:)';
ytmp = ytmp(:)';
%
%  Connectivity.
%
lc1  = length(con(:));
con1 = repmat(con(:)', [1 n2]);
next = (n1+1)*repmat(0:(n2-1), [lc1 1]);
con1 = con1 + next(:)';
con2 = con1 + n1 + 1;
%
n12 = n1*n2;
con1 = reshape(con1, [2 n12]);
con2 = reshape(con2, [2 n12]);
con  = [con1; con2];
%
if (nd == 2)
  crd = [xtmp; ytmp];
  mesh = MeshStructure(crd, con, [], 'ElementType', 'quads:4');
  mesh = ApplyTransformation(mesh, options.Transformation, options.Translation);
  mesh.ndiv = ndiv;
  return
end
%
%  Works for 2D so far ...
%
%  Dimension 3.
%  -----------
%
n3   = ndiv(3);
dz   = 1.0/n3;
zcrd = 0:dz:1;
%
[xtmp, ytmp, ztmp] = ndgrid(xcrd, ycrd, zcrd);
xtmp = xtmp(:)';
ytmp = ytmp(:)';
ztmp = ztmp(:)';
%
%  Connectivity.
%
lc1  = length(con(:));
con1 = repmat(con(:)', [1 n3]);
next = ((n1+1)*(n2+1))*repmat(0:(n3-1), [lc1 1]);
con1 = con1 + next(:)';
con2 = con1 + (n1 + 1)*(n2 + 1);
%
n123 = n1*n2*n3;
con1 = reshape(con1, [4 n123]);
con2 = reshape(con2, [4 n123]);
con  = [con1; con2];
%
crd = [xtmp; ytmp; ztmp];
mesh = MeshStructure(crd, con, [], 'ElementType', 'bricks:8');
mesh = ApplyTransformation(mesh, options.Transformation, options.Translation);
mesh.ndiv = ndiv;
%
%  Now convert to other element type as necessary.
%
if strcmp(options.ElementType, 'tets:4')
  mesh = MeshToSimplex(mesh);
end
%
if strcmp(options.ElementType, 'tets:10')
  mesh = MeshToSimplex(mesh);
  mesh = MeshT4ToT10(mesh);
end
%
return
%
%-------------------- Helper functions.
%
function m = ApplyTransformation(m0, A, x0)
% M - Apply transformation and translation.
%   
ndA  = size(A, 2);
ndx0 = size(x0, 1);
ndm0 = size(m0.crd, 1);
%
if (any([ndA ndx0 ndm0] - ndm0))
  error('dimension mismatch on transformation')
end
%
m = m0;
m.crd = A*m0.crd + repmat(x0, [1 size(m0.crd, 2)]);
%
m.transform.A  = A;
m.transform.x0 = x0;
