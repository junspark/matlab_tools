clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.Image.pname        = '.';
XRDIMAGE.Image.fbase        = 'CeO2_';
XRDIMAGE.Image.fnumber      = 570;
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge1';
XRDIMAGE.Image.corrected    = 0;
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
title('pick three points from a reflection and update script')

p1  = [111 28];
p2  = [722 1420];
p3  = [1838 2028];

x1  = p1(1);
y1  = p1(2);

x2  = p2(1);
y2  = p2(2);

x3  = p3(1);
y3  = p3(2);

LHS = [ ...
    (x1 - x2) (y1 - y2); ...
    (x1 - x3) (y1 - y3); ...
    ];

RHS = 0.5.*[ ...
    (x1*x1 + y1*y1 - x2*x2 - y2*y2); ...
    (x1*x1 + y1*y1 - x3*x3 - y3*y3); ...
    ];

pc = LHS\RHS

xc  = pc(1);
yc  = pc(2);

sqrt((x1 - xc)^2 + (y1 - yc)^2)
sqrt((x2 - xc)^2 + (y2 - yc)^2)
sqrt((x3 - xc)^2 + (y3 - yc)^2)