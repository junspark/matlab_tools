function [] = PlotImage(imdata, imax, imin)
imagesc(imdata)
axis square tight
caxis([imin imax])

% PlotImage(frame_data, max(frame_data(:)), min(frame_data(:)))