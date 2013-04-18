function [delta_tth, delta_eta] = vff_AngularCoverage(r_grid, eta_grid, det_x, det_y, Dsam)
% vff_AngularCoverage - calculates the angular coverage of a panel or a pixel
%
%   USAGE:
%
%   [det_delta_tth, det_delta_eta] = vff_AngularCoverage(r_grid, eta_grid, det_x, det_y, Dsam)
%
%   INPUT:
%
%   r_grid
%       radial grid of the vff detector
%
%   eta_grid
%       azimuthal grid of the vff detector
%
%   det_x
%       x-dimension of the detector
%
%   det_y
%       y-dimension of the detector
%
%   Dsam
%       sample to detector distance
%
%   OUTPUT:
%
%   delta_tth
%       2-theta coverage
%
%   delta_eta
%       eta coverage

%%% SETUP PANEL MESH
crd0    = [ ...
    +det_x/2 +det_x/2 -det_x/2 -det_x/2; ...
    -det_y/2 +det_y/2 +det_y/2 -det_y/2; ...
    Dsam Dsam Dsam Dsam; ...
    ];

%%% CALCULATE
delta_tth   = zeros(length(eta_grid), length(r_grid));
delta_eta   = zeros(length(eta_grid), length(r_grid));
for i = 1:1:length(eta_grid)
    eta = eta_grid(i);
    Rz  = [ ...
        cosd(eta) -sind(eta) 0; ...
        sind(eta) cosd(eta) 0; ...
        0 0 1;
        ];
    for j = 1:1:length(r_grid)
        r   = [r_grid(j) 0 0]';
        
        %%% FIRST PLACE PANEL AT THE CORRECT LOCATION
        crd = Rz*(crd0 + repmat(r, 1, size(crd0,2)));
        
        %%% TTH CALCULATION
        crd_normalized  = NormVecArray(crd);
        crd_normalized  = crd./repmat(crd_normalized,3,1);
        
        tth_min = min(acosd(dot(crd_normalized, repmat([0 0 1]', 1, size(crd0, 2)))));
        tth_max = max(acosd(dot(crd_normalized, repmat([0 0 1]', 1, size(crd0, 2)))));
        
        delta_tth(i,j)  = tth_max - tth_min;
        
        %%% ETA CALCULATION
        crd_det     = crd(1:2,:);
        crd_cen_det = mean(crd_det,2);
        
        crd_det_normalized  = NormVecArray(crd_det);
        crd_det_normalized  = crd_det./repmat(crd_det_normalized, 2, 1);
        
        crd_cen_det_normalized  = NormVecArray(crd_cen_det);
        crd_cen_det_normalized  = crd_cen_det./repmat(crd_cen_det_normalized, 2, 1);
        
        eta_max = mean(acosd(dot(crd_det_normalized, repmat(crd_cen_det_normalized, 1, size(crd0,2)))));
        
        delta_eta(i,j)  = 2*eta_max;
    end
end
