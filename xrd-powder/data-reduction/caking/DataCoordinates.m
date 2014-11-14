function values = DataCoordinates(XY, L, mesh, Ifunc)
% determines the intensity values of the nodal points

crd = mesh.crd;
con = mesh.con;

%square grid
gridNum     = (L-1)*(fix(XY(2))-1)+fix(XY(1));
elemNums    = [gridNum*2-1 gridNum*2];

P   = crd(:, con(:,elemNums(1)));
A1  = [P(1,1)-P(1,3) P(1,2)-P(1,3);P(2,1)-P(2,3) P(2,2)-P(2,3)];
B1  = [XY(1)-P(1,3); XY(2)-P(2,3)];
Y1  = A1\B1;

if all(Y1>=0)
    fele    = elemNums(1);
    fcrd    = [Y1(1) Y1(2) 1-Y1(1)-Y1(2)];
else
    P   = crd(:,con(:,elemNums(2)));
    A2  = [P(1,1)-P(1,3) P(1,2)-P(1,3);P(2,1)-P(2,3) P(2,2)-P(2,3)];
    B2  = [XY(1)-P(1,3); XY(2)-P(2,3)];
    Y2  = A2\B2;
    
    fele    = elemNums(2);
    fcrd    = [Y2(1) Y2(2) 1-Y2(1)-Y2(2)];
end

fcon    = Ifunc(con(:, fele));
values  = dot(fcon, fcrd', 1);