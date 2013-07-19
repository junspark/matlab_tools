clear all
close all
clc

pname   = 'C:\Users\s1iduser\Documents\GitHub\matlab_tools\HEXRD\HEXRDLog_examples';
fname   = 'grainlist_44,45.log';
pfname  = fullfile(pname, fname);

log = parseHEXRDLog(pfname);

% define colors
ncmap   = 256;
cmap    = jet(ncmap);

Completeness    = [log(:).Completeness];
MaxRange    = max(Completeness);
MinRange    = min(Completeness);
dData       = MaxRange - MinRange;
DataColor   = round((ncmap-1)*(Completeness - MinRange)./dData) + 1;


figure,
hold on
for i = 1:1:length(log)
    if ~isnan(Completeness(i))
        xCOM    = log(i).COM(1);
        yCOM    = log(i).COM(2);
        zCOM    = log(i).COM(3);
        
        plot3(xCOM, yCOM, zCOM, ...
            'Marker', 'o', ...
            'MarkerFaceColor', cmap(DataColor(i), :), ...
            'MarkerEdgeColor', cmap(DataColor(i), :), ...
            'MarkerSize', 10);
    end
end
axis square tight
view([46, 44])
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
colorbar vert
caxis([MinRange MaxRange])
grid on

