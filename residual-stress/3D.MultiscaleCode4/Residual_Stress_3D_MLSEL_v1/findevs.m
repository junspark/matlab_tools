function [evals,evecs]=findevs(x,y,z,ndv)

for i=1:ndv
    xc=x(i);
    yc=y(i);
    zc=z(i);
    
    i11=0;i22=0;
    i33=0;
    i12=0;
    i23=0;i13=0;
    
    for j=1:ndv
        % First find the rectangular box that inscribes the data points
        xx=x(j);
        yy=y(j);
        zz=z(j);


        dx=xx-xc;
        dy=yy-yc;
        dz=zz-zc;
        dist=sqrt(dx^2+dy^2+dz^2);
        if not(dist==0)

            % a=3*std(xx);
            % b=3*std(yy);
            % c=3*std(zz);


            i11 = i11+(dx*dx/dist^4);
            i22 = i22+(dy*dy/dist^4);
            i33 = i33+(dz*dz/dist^4);
            i12 = i12+(dx*dy/dist^4);
            i13 = i13+(dx*dz/dist^4);
            i23 = i23+(dz*dy/dist^4);


            tensor=[i11 i12 i13;
                    i12 i22 i23;
                    i13 i23 i33];

        end
        
    end
    

    [evec eval]=eig(inv(tensor));

    a=sqrt(eval(1,1))*2;
    b=sqrt(eval(2,2))*2;
    c=sqrt(eval(3,3))*2;
    
    evals(i,:)=[a b c];
    evecs(i,1:3,1:3)=evec;
    
end



