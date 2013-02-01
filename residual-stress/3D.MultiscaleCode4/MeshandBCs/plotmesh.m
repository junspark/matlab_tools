% Eralp Demir
% 23 Sep 2007
% 10 Nov 2007
% 18 Nov 2007
%
%
%close all
%               
% _______________________
%figure
hold on
%
axis equal
view(3)
% Plot the elements
for i=1:numel
    Node1=np(i,1);
    Node2=np(i,4);
    Node3=np(i,3);
    Node4=np(i,2);
    Node5=np(i,5);
    Node6=np(i,8);
    Node7=np(i,7);
    Node8=np(i,6);
    xx=[x(Node1) x(Node2) x(Node3) x(Node4)...
       x(Node5) x(Node6) x(Node7) x(Node8)];
    yy=[y(Node1) y(Node2) y(Node3) y(Node4)...
       y(Node5) y(Node6) y(Node7) y(Node8)];
    zz=[z(Node1) z(Node2) z(Node3) z(Node4)...
       z(Node5) z(Node6) z(Node7) z(Node8)];
    plot3(xx,yy,zz,'k.')
    line([xx(1:4) xx(1)],[yy(1:4) yy(1)],[zz(1:4) zz(1)],'Color','k');
    line([xx(5:8) xx(5)],[yy(5:8) yy(5)],[zz(5:8) zz(5)],'Color','k');
    line([xx(1) xx(5)],[yy(1) yy(5)],[zz(1) zz(5)],'Color','k');
    line([xx(2) xx(6)],[yy(2) yy(6)],[zz(2) zz(6)],'Color','k');
    line([xx(3) xx(7)],[yy(3) yy(7)],[zz(3) zz(7)],'Color','k');
    line([xx(4) xx(8)],[yy(4) yy(8)],[zz(4) zz(8)],'Color','k');
end

% Plot the surfaces
% for i=1:numels
%     Node1=nps(i,1);
%     Node2=nps(i,2);
%     Node3=nps(i,3);
%     Node4=nps(i,4);
%     xxx=[x(Node1) x(Node2) x(Node3) x(Node4)];
%     yyy=[y(Node1) y(Node2) y(Node3) y(Node4)];
%     zzz=[z(Node1) z(Node2) z(Node3) z(Node4)];
%     line([xxx(1:4) xxx(1)],[yyy(1:4) yyy(1)],[zzz(1:4) zzz(1)],'Color','k');  
% end


% Plot the surfaces
% for i=1:numelsyms
%     Node1=npsyms(i,1);
%     Node2=npsyms(i,2);
%     Node3=npsyms(i,3);
%     Node4=npsyms(i,4);
%     xxx=[x(Node1) x(Node2) x(Node3) x(Node4)];
%     yyy=[y(Node1) y(Node2) y(Node3) y(Node4)];
%     zzz=[z(Node1) z(Node2) z(Node3) z(Node4)];
% %     for kk=1:4
% %         plot3(xxx(kk),yyy(kk),zzz(kk),'b*')
% %         pause
% %     end
% %     plot3(xx,yy,zz,'r*')
% %     pause
% %     line([xxx(1:4) xxx(1)],[yyy(1:4) yyy(1)],[zzz(1:4) zzz(1)],'Color','k');
% end


% 
% % Contour plot of bottom surface of the disk
% figure
% for i=1:55
%     nodesx(i)=x(i);
%     nodesy(i)=y(i);
%     nodesz(i)=gsig(6*i-5);
% end
% % Create a mesh
% DT=delaunay(nodesx,nodesy);
% h=trisurf(DT,nodesx,nodesy,nodesz);
% axis vis3d
% 
% axis off
% l=light('Position',[-50 -15 29]);
% set(gca,'CameraPosition',[208 -50 7687]);
% %lighting phong
% %shading interp
% colorbar EastOutside
%
%
% Plot the stresses in the elements
%figure
%
% % Scaling factor of the results
% scale=0.01;
% % Stress component to be plotted
% %dof=1;
% % sigma-xx
% X=0;Y=0;Z=0;
% for i=1:numnp
%     Node=i;
%     X=[X x(Node)];
%     Y=[Y y(Node)];
%     Z=[Z gsig(6*Node-5)];
% end
% ss=size(X);
% XX=X(2:ss(2));
% YY=Y(2:ss(2));
% ZZ=Z(2:ss(2));
% plot3(XX,YY,ZZ*scale,'r*')
% % sigma-yy
% X=0;Y=0;Z=0;
% for i=1:numnp
%     Node=i;
%     X=[X x(Node)];
%     Y=[Y y(Node)];
%     Z=[Z gsig(6*Node-4)];
% end
% ss=size(X);
% XX=X(2:ss(2));
% YY=Y(2:ss(2));
% ZZ=Z(2:ss(2));
% plot3(XX,YY,ZZ*scale,'b*')
%
% % Make the contour plot
% Xdiv=6;
% Ydiv=numnp/Xdiv/2;
% idx=1;
% idy=1;
% for i=1:numnp
%     if mod(i,Xdiv*2)<=Xdiv & not(mod(i,Xdiv*2)==0)
%         XXX(idx,idy)=XX(i);
%         YYY(idx,idy)=YY(i);
%         ZZZ(idx,idy)=ZZ(i);
%         idx=idx+1;
%     end
%     if mod(i,Xdiv*2)==0
%         idy=idy+1;
%         idx=1;
%     end
% end
% figure
% surf(XXX,YYY,ZZZ)
% colorbar
% figure
% contourf(XXX,YYY,ZZZ)
%
% Plot using the shape functions
% Divide the mapped reference into small segments
%npt=10;
%coord=linspace(0,1,npt);
% Calculate the shape function values at that point and map the coordinates
% using the shape function values
% q3=0;
% counter=0;
% for iel=1:1:numel
%     for  j=1:1:npt
%         q1=coord(j);
%         for  k=1:1:npt
%             q2=coord(k);
%             counter=counter+1;
%             for i=1:1:nnpe
%                 Node=np(iel,i);
%                 xyz_nodes(i,1:3)=[x(Node) y(Node) 0];
%                 sigma_nodes(i,1:3)=[gsig(6*Node-5) gsig(6*Node-4) gsig(6*Node-3)];
%                 sfacplot(i)=sfn3d(q1,q2,q3,i,meltyp);
%             end
%             % Multiply the shape function vector to get the coordinates
%             xyz=sfacplot*xyz_nodes;
%             sigma=sfacplot*sigma_nodes;
%             XXX(counter)=xyz(1);
%             YYY(counter)=xyz(2);
%             ZZZ(counter)=sigma(1);
%         end
%     end
% end
% 
% 
% % Delanuay traingulation for plotting
% % Clear the triangles that have points outside the borders
% DT=delaunay(XXX,YYY);
% dts=size(DT,1);
% DTnew=DT;
% for i=1:dts
%     % Center coordinates of the triangles
%     xc=mean([XXX(DTnew(i,1)) XXX(DTnew(i,2)) XXX(DTnew(i,3))]);
%     yc=mean([YYY(DTnew(i,1)) YYY(DTnew(i,2)) YYY(DTnew(i,3))]);
%     [out] = borders(xc,yc);
%     if out==1
%         DTnew(i,:)=0;
%     end
% end
% index=find(DTnew(:,1));
% DTnew=DTnew(index,:);
% 
% figure
% mesh=trimesh(DTnew,XXX,YYY,ZZZ);
%
% 
% figure
% h=trisurf(DTnew,XXX,YYY,ZZZ);
% axis vis3d
% 
% axis off
% l=light('Position',[-50 -15 29]);
% set(gca,'CameraPosition',[208 -50 7687]);
% %lighting phong
% shading interp
% colorbar EastOutside
% 