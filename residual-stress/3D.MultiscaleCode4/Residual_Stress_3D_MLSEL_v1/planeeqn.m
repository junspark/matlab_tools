% Coordinates
xyz1=[1 1 0.5]';
xyz2=[0.5 1 0]';
xyz3=[1 0.5 0]';

D=[ xyz1(1) xyz1(2) xyz1(3);
    xyz2(1) xyz2(2) xyz2(3);
    xyz3(1) xyz3(2) xyz3(3)];

detD=det(D);

% Coefficients for the plane equation
% x/a + y/b + z/c - 1 = 0
% a, b, c are the intercept values

a=detD/det( [   1 xyz1(2) xyz1(3);
                1 xyz2(2) xyz2(3);
                1 xyz3(2) xyz3(3)]);

b=detD/det( [   xyz1(1) 1 xyz1(3);
                xyz2(1) 1 xyz2(3);
                xyz3(1) 1 xyz3(3)]);
      
c=detD/det( [   xyz1(1) xyz1(2) 1;
                xyz2(1) xyz2(2) 1;
                xyz3(1) xyz3(2) 1]);      