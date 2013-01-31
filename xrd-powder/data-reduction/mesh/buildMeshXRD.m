function mesh=buildMeshXRD(imageData)
L=size(imageData,2);
con=zeros(3,(L-1)*2*(L-1));


for j=1:size(imageData,2)-1
    con(:,1+(L*2-2)*(j-1):(L*2-2)*j)=[[L*(j-1)+1 reshape(repmat((L*(j-1)+2):(L*(j-1)+L-1),2,1),1,(L-2)*2)  L*j];...
        reshape([(L*(j-1)+2):(L*j); (L*j+2):(L*j+L)],1,(L-1)*2);...
        reshape(repmat((L*j+1):(L*j+L-1),2,1),1,(L-1)*2)];
end
x=repmat(1:L,1,L);
y=reshape(repmat(1:L,L,1),1,L^2);
crd=[x;y];
mesh = MeshStructureXRD(crd, con, [], 'ElementType','triangles:3');

load('qr_trid03p06.mat');
mesh.qrule=qr_trid03p06;