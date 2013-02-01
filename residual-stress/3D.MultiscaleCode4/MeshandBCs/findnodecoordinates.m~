function [c,x,y,z]=findnodecoordinates(R,t,alpha)
%
noR=size(R,2);
not=size(t,2);
noalpha=size(alpha,2);
c=0;
for i=1:1:noR
    for j=1:1:noalpha
        for k=1:1:not
            c=c+1;
            x(c)=R(i)*cos(alpha(j)*pi/180);
            y(c)=R(i)*sin(alpha(j)*pi/180);
            z(c)=t(k);
        end
    end
end
