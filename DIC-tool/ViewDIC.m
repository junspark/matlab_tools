clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES
% 1. ASSUMES THAT THE Y DIRECTION (APS CRD) IS HORIZONTAL DIRECTION IN THE DIC
% IMAGE
% 2. EXAMINE CAREFULLY THE COORDINATE SYSTEM(S) BEFORE USING THIS SCRIPT WITH
% IN-SITU LOADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DIC IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname   = 'C:/Users/parkjs/Documents/Projects/Ops/UserSupport/PUP_AFRL_201310/DIC';
fname0  = 'Ti7_00018.tif';          % Reference state
pfname0 = fullfile(pname, fname0);
imdata0 = imread(pfname0);

fname   = 'Ti7_04495.tif';          % Current state
pfname  = fullfile(pname, fname);
imdata  = imread(pfname);

pix2mm  = 0.002;    %%% mm / pixel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ri  = 881;
rf  = 1278;
ci  = 560;
cf  = 2320;
x   = (ci:1:cf) - ci;
x   = x*pix2mm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_box   = [ci cf cf ci ci];
r_box   = [ri ri rf rf ri];

profile0    = imdata0(ri:rf,ci:cf);
profile0    = sum(profile0,1);

profile     = imdata(ri:rf,ci:cf);
profile     = sum(profile,1);

f1  = figure(1);
subplot(1,2,1)
imagesc(imdata0)
colormap(gray)
title(fname0)
axis equal tight off
hold on
line(c_box, r_box)
xlabel('y')
ylabel('x')
colorbar vert
hold off

subplot(1,2,2)
imagesc(imdata)
colormap(gray)
title(fname)
axis equal tight off
hold on
line(c_box, r_box)
xlabel('y')
ylabel('x')
colorbar vert
hold off

f2  = figure(2);
subplot(3,1,1)
imagesc(imdata0(ri:rf,ci:cf))
colormap(gray)
title(['ROI1 - feature in ', fname0])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

subplot(3,1,2)
imagesc(imdata(ri:rf,ci:cf))
colormap(gray)
title(['ROI1 - feature in ', fname])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

subplot(3,1,3)
plot(x, profile0, 'rs-')
hold on
plot(x, profile, 'bo-')
xlabel('y (mm)')
ylabel('summed intensity (arb units)')
