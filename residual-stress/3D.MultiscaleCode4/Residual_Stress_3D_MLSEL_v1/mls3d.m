% Eralp Demir
% Oct. 7, 2010

% This function calculates the stress at the quadrature points using the
% measured stress data at the diffraction volumes using 'Moving Least
% Squares Method'.


function [PHI] = mls3d(xqpt,yqpt,zqpt,dvc_x,dvc_y,dvc_z,rho,iele,nqptv,ndv)


%s=size(dvc_x,2);
% 1D polynomial basis at each node
for i=1:ndv
    % Linear
    p_I(i,:)=[1 dvc_x(i) dvc_y(i) dvc_z(i)];
end


% At the quad points in the domain (x)
%csig=zeros(6,nqptv,numel);
%for iele=1:numel
    for iquad=1:1:nqptv
        x=xqpt(iele,iquad);
        y=yqpt(iele,iquad);
        z=zqpt(iele,iquad);
        A=zeros(4,4);
        
        for i=1:ndv
            r=sqrt((x-dvc_x(i))^2+(y-dvc_y(i))^2+(z-dvc_z(i))^2);
            w=mls3dweight(r,rho);
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
             disp('Rank deficiency!!!');
             disp('Problem occurs in MLS method!!!');
             disp('************************************************************');
             % phi is cleared to generate an error and stop running!!!
             clear phi
         else
            phi=[1 x y z]*inv(A)*B;        
         end
         
        PHI(iquad,1:ndv)=phi;
    end
%end



