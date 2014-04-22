function mesh = BuildMeshDetector(Lx, Ly)
% BuildMeshDetector2 - Build mesh crd and con for an XRD image
%   
%   mesh = MeshStructure2(Lx, Ly)
%
%   Lx is the number of pixels in the image in the horizontal direction.
%   Ly is the number of pixels in the image in the vertical direction.
%
%   mesh is a structure with appropriate fields
%
if nargin == 1
    disp('square image')
    Ly = Lx;
end

numel   = (Lx - 1) * (Ly - 1) * 2;
con     = zeros(3, numel);

for i = 1:1:(Ly - 1)
    idx1    = 1 + 2*(Lx - 1) * ( i - 1 );
    idx2    = 2*(Lx - 1) * i;
    
    xL  = Lx*(i - 1)+1;
    xR  = Lx*i;
    xLR = (xL + 1):1:(xR - 1);
    xLR = repmat(xLR, 2, 1);
    
    con(1, idx1: idx2)  = [xL xLR(:)' xR];
    
    xL  = Lx*(i - 1) + 2 : 1 : Lx * i;
    xR  = Lx*i + 2 : 1 : Lx * (i + 1);
    xLR = [xL; xR];
    
    con(2, idx1: idx2)  = xLR(:)';
    
    xL  = Lx * i + 1;
    xR  = Lx * (i + 1) - 1;
    xLR = xL:xR;
    xLR = repmat(xLR, 2, 1);
    con(3, idx1: idx2)  = xLR(:)';
    
end
x   = repmat(1:Lx, 1, Ly);
y   = repmat(1:Ly, Lx, 1);
y   = y(:)';

crd = [x;y];

load('qr_trid03p06.mat');
mesh = MeshStructureXRD(crd, con, numel, qr_trid03p06);