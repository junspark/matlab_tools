function Akeep = BuildLsdfSpfMatrix(hkl, mesh, sym, pts, div, odf, w, deta, dome)
% BuildLsdfSpfMatrix - Build ODF/PF matrix in pieces.
%   
%   USAGE:
%
%   Akeep = BuildLsdfSpfMatrix(hkl, mesh, sym, pts, div, odf, w)
%
%   INPUT:
%
%   hkl, mesh, sym, pts, div, odf, w
%          are the same as in `WeightedOdfPfMatrix'
%
%   OUTPUT:
%
%   Akeep is output from `WeightedOdfPfMatrix' with projection operator
%
%   NOTES:
%
%   *  See WeightedOdfPfMatrix for further documentation.
%

%%%
pts     = FiberTube(pts, deta, dome);
npts    = size(pts, 2); 

%
nnode   = mesh.numind;

Akeep   = zeros(npts, nnode*6);
for ii = 1:1:npts
    Aopm    = WeightedOdfPfMatrix(hkl, mesh, sym, pts(:,ii), div, odf);
    Gamma   = ProjectionVector(pts(:,ii));
    
    G   = zeros(mesh.numind, 6*mesh.numind);
    for j = 1:1:mesh.numind
        ci  = 6*(j-1) + 1;
        cf  = 6*j;
        G(j,ci:cf)  = Gamma;
    end
    A   = Aopm*G;
    Akeep(ii,:) = A;
end

Akeep   = sparse(Akeep);
Akeep   = w*Akeep;