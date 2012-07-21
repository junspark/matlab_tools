function Akeep = BuildLsdfSpfMatrix(hkl, mesh, sym, pts, div, odf, w)
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
[TH,PHI,~] = cart2sph(pts(1,1),pts(2,1),pts(3,1));

V1  = deg2rad(2.5);
V2  = deg2rad(2.5)*cosd(45);

THpts   = [TH+V1; TH+V2; TH; TH-V2; TH-V1; TH-V2; TH; TH+V2];
PHIpts  = [PHI; PHI+V2; PHI+V1; PHI+V2; PHI; PHI-V2; PHI-V1; PHI-V2];
Rpts    = ones(size(THpts));

[x, y, z]   = sph2cart(THpts,PHIpts,Rpts);

pts = [pts [x y z]'];
pts = UnitVector(pts);

%
npts    = size(pts, 2); 
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