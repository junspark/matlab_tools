mls_sig=reshape(gsig,6,numnp);

% % Plot the 1st component of stress at all nodes
% figure
% hold on
% plot3(x,y,mls_sig(1,:),'r*');
% plot3(x,y,as_sig(1,:),'b+');
% %
for i=1:size(x,2)
    r(i)=sqrt(x(i)^2+y(i)^2);
end
figure
hold on
plot(r,mls_sig(3,:),'b*');
%plot(r,as_sig(1,:),'k+');