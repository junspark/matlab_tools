clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.Image.pname        = 'C:\Users\parkjs\Documents\GitHub\matlab_tools_examples\xrd-powder-data-reduction-example\APS\hydra_example_meimei_aug14\';
XRDIMAGE.Image.fbase        = 'ceria_';
XRDIMAGE.Image.fnumber      = 15;
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge1.sum';
XRDIMAGE.Image.corrected    = 1;
XRDIMAGE.Image.IsHydra      = 1;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% DARK FILES ONLY USED IF THE IMAGES ARE UNCORRECTED
XRDIMAGE.DarkField.pname    = '.';
XRDIMAGE.DarkField.fbase    = 'dark_';
XRDIMAGE.DarkField.fnumber  = 571;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 5;
XRDIMAGE.DarkField.fext     = 'ge1';

XRDIMAGE.Calib.pname        = '.';
XRDIMAGE.Calib.fbase        = 'CeO2_';
XRDIMAGE.Calib.fnumber      = 570;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD XRD IMAGES
%%% BACKGROUND
if XRDIMAGE.Image.corrected
    disp('###########################')
    fprintf('images are already corrected for background.\n');
    disp('###########################')
else
    disp('###########################')
    fprintf('loading background file for dark.\n');
    disp('###########################')
    pfname  = GenerateGEpfname(XRDIMAGE.DarkField);
    bg      = NreadGE(pfname{1,1}, 1);
end

pfname  = GenerateGEpfname(XRDIMAGE.Image);
i   = 1;
disp('###########################')
disp(sprintf('Looking at %s', pfname{i,1}))
disp('###########################')
if XRDIMAGE.Image.corrected
    imgi    = ReadSUM(pfname{i,1});
else
    imgi    = bg.*0;
    for j = 1:1:XRDIMAGE.Image.numframe
        imgj    = NreadGE(pfname{i,1}, j);
        imgi    = imgi + imgj;
    end
    imgi    = imgi - bg.*XRDIMAGE.Image.numframe;
end

figure(1)
hold off
imagesc(rot90(imgi,1))
caxis([-10 3000])
axis equal
colorbar vert
hold on
xlabel('X_L (pixels)')
ylabel('Y_L (pixels)')
title('pick three points from a reflection')

[x, y] = ginput(3);

x1  = x(1);
y1  = y(1);

x2  = x(2);
y2  = y(2);

x3  = x(3);
y3  = y(3);

LHS = [ ...
    (x1 - x2) (y1 - y2); ...
    (x1 - x3) (y1 - y3); ...
    ];

RHS = 0.5.*[ ...
    (x1*x1 + y1*y1 - x2*x2 - y2*y2); ...
    (x1*x1 + y1*y1 - x3*x3 - y3*y3); ...
    ];

pc = LHS\RHS;

xc  = pc(1);
yc  = pc(2);

sqrt((x1 - xc)^2 + (y1 - yc)^2)
sqrt((x2 - xc)^2 + (y2 - yc)^2)
sqrt((x3 - xc)^2 + (y3 - yc)^2)

disp('copy the following lines to appropriate location in the subsequent analysis script')
disp(sprintf('XRDIMAGE.CakePrms.origin(1) = %f;         %% apparent X center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1', pc(1)))
disp(sprintf('XRDIMAGE.CakePrms.origin(2) = %f;            %% apparent Y center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1', pc(2)))