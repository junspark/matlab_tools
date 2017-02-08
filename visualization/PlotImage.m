function PlotImage(image_data, max_range, min_range)
figure,
clf
%%% CHANGE RANGE AS DESIRED
imagesc(image_data, [min_range max_range])
colorbar vert
axis square tight