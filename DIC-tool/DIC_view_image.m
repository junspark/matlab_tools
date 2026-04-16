clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIC_view_image.m
% Derived from ViewDIC.m
%
% NOTES
% 1. ASSUMES THAT THE Y DIRECTION (APS CRD) IS HORIZONTAL DIRECTION IN THE DIC
% IMAGE
% 2. EXAMINE CAREFULLY THE COORDINATE SYSTEM(S) BEFORE USING THIS SCRIPT WITH
% IN-SITU LOADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXPERIMENT SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_info_file = 'dic_exp_info.md';     % path to experiment config file

[pname, froot, fext, ndigits, fini, ffin, finc, ...
    delta_h, delta_v, Hctr, Vctr, ws_h, ws_v, pix2mm, gauge_length] = dic_read_exp_setup(exp_info_file);

ci  = Hctr - ws_h * delta_h;   % ROI start column (pixels)
cf  = Hctr + ws_h * delta_h;   % ROI end column (pixels)
ri  = Vctr - ws_v * delta_v;   % ROI start row (pixels)
rf  = Vctr + ws_v * delta_v;   % ROI end row (pixels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DIC IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname0  = sprintf('%s_%s.%s', froot, sprintf(['%0', num2str(ndigits), 'd'], fini), fext);
fname   = sprintf('%s_%s.%s', froot, sprintf(['%0', num2str(ndigits), 'd'], ffin), fext);

pfname0 = fullfile(pname, fname0);
pfname  = fullfile(pname, fname);

imdata0 = imread(pfname0);
imdata  = imread(pfname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x   = (ci:1:cf) - ci;
x   = x * pix2mm;

c_box   = [ci cf cf ci ci];
r_box   = [ri ri rf rf ri];

profile0    = imdata0(ri:rf, ci:cf);
profile0    = sum(profile0, 1);

profile     = imdata(ri:rf, ci:cf);
profile     = sum(profile, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
imagesc(imdata0(ri:rf, ci:cf))
colormap(gray)
title(['ROI1 - feature in ', fname0])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

subplot(3,1,2)
imagesc(imdata(ri:rf, ci:cf))
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
