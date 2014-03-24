function Ilist = PolarIntensityXRD(R,V,delTH,qrule)
% Rlist=zeros(size(V,2)-1,1);
Ilist=zeros(size(V,2)-1,1);

for jj=1:size(V,1)-1
    VL=zeros(size(V,2)-1,1);
    VR=zeros(size(V,2)-1,1);
    AL=zeros(size(V,2)-1,1);
    AR=zeros(size(V,2)-1,1);
%     Rval=zeros(size(V,2)-1,1);
    
    for ii=1:size(V,2)-1
        vrwL=zeros(1,size(qrule.pts,2));
        ArwL=zeros(1,size(qrule.pts,2));
        vrwR=zeros(1,size(qrule.pts,2));
        ArwR=zeros(1,size(qrule.pts,2));
        for kk=1:size(qrule.pts,2)
            v=V(jj,ii)*qrule.pts(1,kk)+V(jj,ii+1)*qrule.pts(2,kk)+V(jj+1,ii+1)*qrule.pts(3,kk);
            r=R(jj,ii)*qrule.pts(1,kk)+R(jj,ii+1)*qrule.pts(2,kk)+R(jj+1,ii+1)*qrule.pts(3,kk);
            Jac=abs(det([R(jj,ii+1)-R(jj,ii) R(jj+1,ii+1)-R(jj,ii); delTH delTH]));
            vrwL(kk)=v*r*qrule.wts(kk)*Jac;
            ArwL(kk)=r*qrule.wts(kk)*Jac;
            
            v=V(jj,ii)*qrule.pts(1,kk)+V(jj+1,ii+1)*qrule.pts(2,kk)+V(jj+1,ii)*qrule.pts(3,kk);
            r=R(jj,ii)*qrule.pts(1,kk)+R(jj+1,ii+1)*qrule.pts(2,kk)+R(jj+1,ii)*qrule.pts(3,kk);
            Jac=abs(det([R(jj+1,ii+1)-R(jj,ii) R(jj+1,ii)-R(jj,ii); delTH 0]));
            vrwR(kk)=v*r*qrule.wts(kk)*Jac;
            ArwR(kk)=r*qrule.wts(kk)*Jac;
        end
        VL(ii)=sum(vrwL);
        AL(ii)=sum(ArwL);
        VR(ii)=sum(vrwR);
        AR(ii)=sum(ArwR);
%         Rval(ii)=mean([R(jj,ii) R(jj+1,ii)]);
        
    end
    
    Ilist(jj)=(sum(VL)+sum(VR))/(sum(AL)+sum(AR));
%     Rlist(jj)=mean(Rval);
end