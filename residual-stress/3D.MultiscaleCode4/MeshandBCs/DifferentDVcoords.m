load tR

alpha=[0 22.5 45 67.5 90];
cc=0;
for i=1:5
    ang=pi*alpha(i)/180;
    for j=1:size(R,1)
        x=R(j)*cos(ang);
        y=R(j)*sin(ang);
        for k=1:size(t,1)
            cc=cc+1;
            z=t(k);
            dvc_x(cc)=x;
            dvc_y(cc)=y;
            dvc_z(cc)=z;
        end
   end
end

plot3(dvc_x,dvc_y,dvc_z,'r*')