function [xqpt,yqpt,zqpt] = qptloc(x,y,z,np,numel,nnpe,sfac)


% This part should be performed before calculations 
% The quadrature point coordinates are calculated using the shape functions
% given by the variable "sfac"


for iele=1:1:numel
    % Calculate the nodal coordinates
    for j = 1:1:nnpe
        j1 = np(iele,j);
        xloc(j) = x(j1);
        yloc(j) = y(j1);
        zloc(j) = z(j1);
    end

    % Calculate the quadrature point coordinates
    xqpt(iele,:)=sfac'*xloc';
    yqpt(iele,:)=sfac'*yloc';
    zqpt(iele,:)=sfac'*zloc';
    
end