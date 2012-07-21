function rmat = RMatOfRoe(roe)
% RMATOFROE - Rotation matrix from Roe (Matthies) Euler angles.
%   
%   rmat = RMatOfRoe(roe)
%
%   roe is 3 x n, the array of Roe parameters (in radians)
%   rmat is 3 x 3 x n, the corresponding rotation matrices
%
% ref: Matthies, S., Vinel, G. W., Helming, K., 
%       "Standard Distributions in Texture Analysis". Vol. 1
%       Akademie-Verlag, Berlin, 1987.  pg. 7

n    = size(roe, 2);
croe = cos(roe);
sroe = sin(roe);

rmat = [...
        croe(1, :).*croe(2, :).*croe(3, :) - sroe(1, :).*sroe(3, :);...
        -1.0*croe(1, :).*croe(2, :).*sroe(3, :) - sroe(1, :).*croe(3, :);...
        croe(1, :).*sroe(2, :);...
        
        sroe(1, :).*croe(2, :).*croe(3, :) + croe(1, :).*sroe(3, :);...
        -1.0*sroe(1, :).*croe(2, :).*sroe(3, :) + croe(1, :).*croe(3, :);...
        sroe(1, :).*sroe(2, :);...
        
        -1.0*sroe(2, :).*croe(3, :);...
        sroe(2, :).*sroe(3, :);...
        croe(2, :);...
    ];

rmat = reshape(rmat, [3 3 n]);

%% y = [0 1 0]';
%% z = [0 0 1]';
%% 
%% rm1 = expm(SkewMatrixOfVector(roe(1)*z));
%% rm2 = expm(SkewMatrixOfVector(roe(2)*y));
%% rm3 = expm(SkewMatrixOfVector(roe(3)*z));
%% 
%% rmat = rm3*rm2*rm1;
