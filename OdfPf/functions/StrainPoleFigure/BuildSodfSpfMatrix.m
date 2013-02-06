function Akeep = BuildSodfSpfMatrix(hkl, mesh, sym, pts, div, odf, w, S, deta, dome)
% BuildSodfSpfMatrix - Build ODF/PF matrix in pieces.
%   
%   USAGE:
%
%   Akeep = BuildSodfSpfMatrix(hkl, mesh, sym, pts, div, odf, w)
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
    
    %%% THIS CAN MOVE OUT OF ii-LOOP
    for jj = 1:1:nnode
        ci  = 1 + 6*(jj - 1);
        cf  = 6*jj;
        R   = RMatOfQuat(QuatOfRod(mesh.crd(:,jj)));
        T   = VectorizedCOBMatrix(R);
        S_STAR  = T*S*T';
        
        Gamma2(jj, ci:cf)   = Gamma*S_STAR;
    end
    A   = Aopm*Gamma2;
    Akeep(ii,:) = A;
end

Akeep   = sparse(Akeep);
Akeep   = w*Akeep;