clear
clc

ang=linspace(0,pi/2,6);

ro=15;

for i=1:6
    xcoord(i)=ro*cos(ang(i));
    ycoord(i)=ro*sin(ang(i));
end
    
for i=1:5
    xmed(i)=mean([xcoord(i) xcoord(i+1)]);
    ymed(i)=mean([ycoord(i) ycoord(i+1)]);
    angmed(i)=mean([ang(i) ang(i+1)]);
    xangmed(i)=ro*cos(angmed(i));
    yangmed(i)=ro*sin(angmed(i));
end

deltax=xangmed-xmed;
deltay=yangmed-ymed;

