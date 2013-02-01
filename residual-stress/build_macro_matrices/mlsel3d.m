% Eralp Demir
% Oct. 7, 2010
%
% This function calculates the stress at the quadrature points using the
% measured stress data at the diffraction volumes using 'Moving Least
% Squares Method'.
%
%
%
function [PHI] = mlsel3d(xqpt,yqpt,zqpt,dvc_x,dvc_y,dvc_z,evals,evecs,scale,iele,nqptv,ndv)
%
%
% a,b,c are the focal points of the ellipse
a=evals(:,1)*scale;
b=evals(:,2)*scale;
c=evals(:,3)*scale;
% evec stores the orientation of the ellipse
% dvc_x, dvc_y, dvc_z are the center points of the ellipse
%
% You must first transform the point under concern into the coordinate
% reference of the ellipsoid and place the coordinate into the ellips
% equation to find 'rho'
%
%
%s=size(dvc_x,2);
% 1D polynomial basis at each node
for i=1:ndv
    % Linear
    p_I(i,:)=[1 dvc_x(i) dvc_y(i) dvc_z(i)];
end
%
%
% At the quad points in the domain (x)
%
%csig=zeros(6,nqptv,numel);
%for iele=1:numel
for iquad=1:1:nqptv
    % Quad points in the global coordinate reference
    x=xqpt(iele,iquad);
    y=yqpt(iele,iquad);
    z=zqpt(iele,iquad);
    A=zeros(4,4);
    %
    for i=1:ndv
        % Vector of the point under concern in the global reference
        vec=[x-dvc_x(i); y-dvc_y(i); z-dvc_z(i)];
        % Orientation of the ellipse
        orient(:,:)=evecs(i,:,:);
        w=mlsel3dweight(vec,orient,a(i),b(i),c(i));
        if i==1
            B=w*p_I(i,:)';
        else
            B=[B w*p_I(i,:)'];
        end
        A=A+(w*p_I(i,:)'*p_I(i,:));
    end
     %
     %
     % Check of rank: will give an error message and calculations will
     % stop HERE!!!
     if rank(A)<4
         disp('************************************************************');
         disp(['element no: ' num2str(iele)]);
         disp(['quad. pt. no: ' num2str(iquad)]);
         disp('Rank deficiency!!!');
         disp('Problem occurs in MLS method!!!');
         disp('************************************************************');
         clear phi
     else
        phi=[1 x y z]*inv(A)*B;
     end
    %
    PHI(iquad,1:ndv)=phi;
end
%end
%
%
%
%