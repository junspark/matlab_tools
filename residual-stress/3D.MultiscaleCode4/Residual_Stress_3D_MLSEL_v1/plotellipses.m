function plotellipses(evals,evecs,i,x,y,z,scale)
%figure 
%hold on
%view([0 0])
%
ndiv=20;
theta=linspace(0,2*pi,ndiv);
phi=linspace(0,pi,ndiv);
%
%for i=1:ndv
   xc=x(i);
   yc=y(i);
   zc=z(i);
    
   a=evals(i,1);
   b=evals(i,2);
   c=evals(i,3);
%
   evec(:,:)=evecs(i,:,:);
%    
    % Plot the ellipsoid
    cc=0;
    for th=1:ndiv
        for ph=1:ndiv
            cc=cc+1;
            xel=scale*a*cos(theta(th))*sin(phi(ph));
            yel=scale*b*sin(theta(th))*sin(phi(ph));
            zel=scale*c*cos(phi(ph));
            xyzel_tr=(evec*[xel;yel;zel])+[xc;yc;zc];
            %xyel_tr=(evec*[xel;yel])+[xc;yc];
            % Plot the points
            %plot3(xyzel_tr(1),xyzel_tr(2),xyzel_tr(3),'r*');
            %plot(xyel_tr(1),xyel_tr(2),'r*');
            xs(cc)=xyzel_tr(1);
            ys(cc)=xyzel_tr(2);
            zs(cc)=xyzel_tr(3);
        end
    end
    plot3(xs,ys,zs,'b');
%
%end
%
%
axis equal