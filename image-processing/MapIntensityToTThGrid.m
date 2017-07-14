function [tth_grid, d_grid, intensity_in_tth_grid] = MapIntensityToTThGrid(XRDIMAGE, polimg)
% MapIntensityToTThGrid - Apply geometric model and detector distortion to
% generate intensity vs. tth data in a standard tth grid

% save('MapIntensityToTThGrid.mat')
% return
% clear all
% close all
% clc
% load('MapIntensityToTThGrid.mat')

tth_min = max(polimg.mapped_tth_for_intensity(:,1));
tth_max = min(polimg.mapped_tth_for_intensity(:,end));
tth_grid    = linspace(tth_min, tth_max, XRDIMAGE.CakePrms.bins(2));
d_grid      = (XRDIMAGE.Instr.wavelength/2)./sind(tth_grid./2);

% intensity_in_tth_grid   = zeros(XRDIMAGE.CakePrms.bins(1), XRDIMAGE.CakePrms.bins(2));
% for i = 1:1:XRDIMAGE.CakePrms.bins(1)
%     for j = 1:1:XRDIMAGE.CakePrms.bins(2)
%         idx1    = find(tth_grid(j) >= polimg.mapped_tth_for_intensity(i,:));
%         idx2    = find(tth_grid(j) <= polimg.mapped_tth_for_intensity(i,:));
%         
%         idx1    = idx1(1);
%         idx2    = idx2(1);
%         
%         dL  = tth_grid(j) - polimg.mapped_tth_for_intensity(i,idx1);
%         dR  = polimg.mapped_tth_for_intensity(i,idx2) - tth_grid(j);
%         
%         intensity_in_tth_grid(i,j)  = (dR*polimg.intensity(i,idx1) + dL*polimg.intensity(i,idx2))/(dL + dR);
%     end
% end
% keyboard

for i = 1:1:XRDIMAGE.CakePrms.bins(1)
    intensity_in_tth_grid2(i,:) = interp1(polimg.mapped_tth_for_intensity(i,:), polimg.intensity(i,:), tth_grid, 'nearest');
end
intensity_in_tth_grid   = intensity_in_tth_grid2;


% imagesc(intensity_in_tth_grid'-intensity_in_tth_grid2')
