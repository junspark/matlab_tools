function [tth_grid, intensity_in_tth_grid] = MapIntensityToTThGrid(XRDIMAGE, polimg)
tth_min = max(polimg.mapped_tth_for_intensity(:,1));
tth_max = min(polimg.mapped_tth_for_intensity(:,end));
tth_grid    = linspace(tth_min, tth_max, XRDIMAGE.CakePrms.bins(2));

intensity_in_tth_grid   = zeros(XRDIMAGE.CakePrms.bins(1), XRDIMAGE.CakePrms.bins(2));
for i = 1:1:XRDIMAGE.CakePrms.bins(1)
    for j = 1:1:XRDIMAGE.CakePrms.bins(2)
        idx1    = find(tth_grid(j) >= polimg.mapped_tth_for_intensity(i,:));
        idx2    = find(tth_grid(j) <= polimg.mapped_tth_for_intensity(i,:));
        
        idx1    = idx1(1);
        idx2    = idx2(1);
        
        dL  = tth_grid(j) - polimg.mapped_tth_for_intensity(i,idx1);
        dR  = polimg.mapped_tth_for_intensity(i,idx2) - tth_grid(j);
        
        intensity_in_tth_grid(i,j)  = (dR*polimg.intensity(i,idx1) + dL*polimg.intensity(i,idx2))/(dL + dR);
    end
end

