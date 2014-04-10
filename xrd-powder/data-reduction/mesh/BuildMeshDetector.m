function mesh = BuildMeshXRD(L)
% BuildMeshXRD -Build mesh crd and con for an XRD image
%   
%   mesh = MeshStructure(L)
%
%   L is the number of pixels in the image in the horizontal direction.
%   This assumes that the image is square and the pixel size in the
%   horizontal and vertical directions are equal
%
%   mesh is a structure with appropriate fields
%
%
numel   = (L-1)*(L-1)*2;
con     = zeros(3,numel);

for j=1:L-1
    con(:,1+(L*2-2)*(j-1):(L*2-2)*j)=[...
        [L*(j-1)+1 reshape(repmat((L*(j-1)+2):(L*(j-1)+L-1),2,1),1,(L-2)*2)  L*j];...
        reshape([(L*(j-1)+2):(L*j); (L*j+2):(L*j+L)],1,(L-1)*2);...
        reshape(repmat((L*j+1):(L*j+L-1),2,1),1,(L-1)*2)];
end
x=repmat(1:L,1,L);
y=reshape(repmat(1:L,L,1),1,L^2);
crd=[x;y];

load('qr_trid03p06.mat');
mesh = MeshStructureXRD(crd, con, numel, qr_trid03p06);