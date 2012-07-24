function [Sagg] = AggregateElasticityMatrix(frmesh, odf, Sxstal)
odf     = odf(:);
nnode   = frmesh.numind;

S   = zeros(6,6,nnode);
for i = 1:1:nnode
    R   = RMatOfQuat(QuatOfRod(frmesh.crd(:,i)));
    T   = VectorizedCOBMatrix(R);
    S(:,:,i)    = T*Sxstal*T';
end

Sagg    = zeros(6,6);
for i = 1:1:6
    for j = 1:1:6
        Sij = S(i,j,:);
        Sij = Sij(:);
        
        Sagg(i,j)   = odf'*frmesh.l2ip*Sij/sum(frmesh.l2ip*odf);
    end
end