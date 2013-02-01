clear all
close all
clc

load inputdata.A1.mat
Alpha   = [0 30 60 90];

% load inputdata.mat
PFNAME_GRID = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_GRID);

nx  = length(inputdata.x_grid);
ny  = length(inputdata.y_grid);
nz  = length(inputdata.z_grid);

ct  = 1;
for i = 1:1:nx
    for j = 1:1:ny
        for ii = 1:1:length(Alpha)
            for k = 1:1:nz
                x   = inputdata.x_grid(i);
                y   = -inputdata.y_grid(j)+6.35;
                z   = inputdata.z_grid(k);
                xyz = [x; y; z];
                
                a   = Alpha(ii);
                R   = [
                    cosd(a) -sind(a) 0;
                    sind(a)  cosd(a) 0;
                    0 0 1;
                    ];
                xyz = R*xyz;
                
                dvc_x(1,ct) = xyz(1);
                dvc_y(1,ct) = xyz(2);
                dvc_z(1,ct) = xyz(3);
                ct  = ct + 1;
            end
        end
    end
end
ndv     = ct - 1;
disp('generating DV grid file')
save(PFNAME_GRID, 'dvc_x', 'dvc_y', 'dvc_z', 'ndv')

figure(2)
for i = 1:1:ndv
    plot3(dvc_x(i), dvc_y(i), dvc_z(i), 'k.')
    axis equal tight
    grid on
    hold on
    pause(0.5)
end