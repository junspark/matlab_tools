function [ref, con] = SubdivideSimplex(dim, n)
% SubdivideSimplex - Regular subdivision of a simplex.
%
%   USAGE:
%
%   ref = SubdivideSimplex(dim, n)
%   [ref, con] = SubdivideSimplex(dim, n)
%
%   INPUT:
%
%   dim is a positive integer, 
%       the dimension of the simplex
%   n   is a positive integer, 
%       the number of subdivisions in each direction
%
%   OUTPUT:
%
%   ref is d x m, 
%       the list of barycentric coordinates of the points 
%       in the subdivision
%   con is d x e, (integer)
%       the connectivity of the subdivision
%
if (dim == 1)
  ref = subdiv_1(n);
  con = connect_1(n);
elseif (dim == 2)
  ref = subdiv_2(n);
  con = connect_2(n);
elseif (dim == 3)
  ref = subdiv_3(n);
  con = connect_3(n, ref);
else
  error(['argument "dim" must be 1, 2, or 3' ...
	 sprintf('\n') ...
	 '   (current value:  dim = ' num2str(dim) ')'])
end
%
%  Fix Jacobians in 3D by reversing first two rows.
%
if (dim == 3)
  mesh = MeshStructure(ref, con);
  jac  = Jacobian(mesh);
  fix  = (jac < 0);
  tmp  = con(1, fix);
  con(1, fix) = con(2, fix);
  con(2, fix) = tmp;
end
%  Return barycentric coordinates.
%
ref = [ref; 1 - sum(ref, 1)];
%
%--------------------Nodal Points by Dimension--------------------------
%                    1D-------------------------------------------------
%
function ref = subdiv_1(n)
%   
ref =(0:n)/n; % order by increasing x1
%
%                    2D-------------------------------------------------
%
function ref = subdiv_2(n)
%   
dim    = 2;
subfun = @subdiv_1;
%
rp = [];
x1 = (0:(n-1))./n; % order by increasing x1
for i=1:n
  ndv  = n + 1 - i;
  sec  = feval(subfun, ndv);
  nsec = size(sec, 2);
  x1i  = x1(i);
  rpi  = [repmat(x1i, [1 nsec]); (1-x1i)*sec];
  rp   = [rp rpi];
end
rp = [rp eye(dim, 1)];
%
ref = rp;
%
%                    3D-------------------------------------------------
%
function ref = subdiv_3(n)
%   
dim    = 3;
subfun = @subdiv_2;
%
rp = [];
x1 = (0:(n-1))./n; % order by increasing x1
for i=1:n
  ndv  = n + 1 - i;
  sec  = feval(subfun, ndv);
  nsec = size(sec, 2);
  x1i  = x1(i);
  rpi  = [repmat(x1i, [1 nsec]); (1-x1i)*sec];
  rp   = [rp rpi];
end
rp = [rp eye(dim, 1)];
%
ref = rp;
%
%--------------------Connectivity by Dimension--------------------------
%                    1D-------------------------------------------------
%
function con = connect_1(n)
%
con = [1:n; 2:(n+1)];
%
%                    2D-------------------------------------------------
%
function con = connect_2(n)
%
np_0 = 1;
con  = [];
for i=1:n
  nsec = n + 2 - i;
  %
  pts3 = np_0:(np_0 + n - i);
  pts1 = pts3 + nsec;
  pts2 = pts3 + 1;
  con1 = [pts1; pts2; pts3];
  %
  con = [con con1];
  %
  pts1 = (np_0+1):(np_0 + n - i);
  pts2 = pts1 + nsec - 1;
  pts3 = pts2 + 1;
  con1 = [pts1; pts2; pts3];
  %
  con = [con con1];
  %
  np_0 = np_0 + nsec;
  %
end
%
%                    3D-------------------------------------------------
%
function con = connect_3(n, r)
%
%  Uses Delaunay.
%
rsum = sum(r);
rsca = 1 + n*rsum;
rsca = repmat(rsca, [3 1]);
rsca(1,:) = 0.98*rsca(1,:);
rsca(2,:) = 0.99*rsca(2,:);
r = r.* rsca;
c = delaunayn(r');
con = c';

