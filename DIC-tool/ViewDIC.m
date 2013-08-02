clear all
close all
clc

% pname   = '/net/s1dserv/export/s1-idb/park_jul2013/DIC4Jun/DIC';
pname   = 'W:\park_jul2013\DIC4Jun\DIC';
fname0  = 'DIC_00045.tif'; % Initial state
% fname0  = 'DIC_00506.tif'; % Final state in the loadframe
pfname0 = fullfile(pname, fname0);
imdata0 = imread(pfname0);

fname   = 'DIC_00101.tif';      %%%%%%%%%%%%%%%%
pfname  = fullfile(pname, fname);
imdata  = imread(pfname);

ri  = 593;
rf  = 709;
ci  = 150;
cf  = 870;
x   = (ci:1:cf) - ci;
x   = x*0.0019;

profile0    = imdata0(ri:rf,ci:cf);
profile0    = sum(profile0,1);

profile     = imdata(ri:rf,ci:cf);
profile     = sum(profile,1);

% f1  = figure(1);
% % set(f1, 'position', [680 678 560 420])
% subplot(2,2,1)
% imagesc(imdata0)
% colormap(gray)
% title(fname0)
% axis equal tight off
% xlabel('y')
% ylabel('x')
% 
% subplot(2,2,2)
% imagesc(imdata)
% colormap(gray)
% title(fname)
% axis equal tight off
% xlabel('y')
% ylabel('x')
% 
% subplot(2,2,3)
% imagesc(imdata0 - imdata)
% colormap(gray)
% title('difference image')
% axis equal tight off
% xlabel('y')
% ylabel('x')
% 
% subplot(2,2,4)
% plot(x, profile0, 'rs-')
% hold on
% plot(x, profile, 'bo-')
% xlabel('y (mm)')
% ylabel('edge (mm)')

ri  = 424;
rf  = 564;
ci  = 2252;
cf  = 2419;
f2  = figure(2);
set(f2, 'position', [1249 178 560 420])
subplot(2,1,1)
imagesc(imdata0(ri:rf,ci:cf))
colormap(gray)
title(['BOTTOM ROI - feature in ', fname0])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

subplot(2,1,2)
imagesc(imdata(ri:rf,ci:cf))
colormap(gray)
title(['BOTTOM ROI - feature in ', fname])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

% ri  = 675;
% rf  = 723;
% ci  = 1369;
% cf  = 1555;
% f3  = figure(3);
% set(f3, 'position', [681 178 560 420])
% subplot(2,1,1)
% imagesc(imdata0(ri:rf,ci:cf))
% colormap(gray)
% title(['BEAM ROI - feature in ', fname0])
% axis equal tight off
% xlabel('y (pixels)')
% ylabel('x (pixels)')
% 
% subplot(2,1,2)
% imagesc(imdata(ri:rf,ci:cf))
% colormap(gray)
% title(['BEAM ROI - feature in ', fname])
% axis equal tight off
% xlabel('y (pixels)')
% ylabel('x (pixels)')
% 
% ri  = 541;
% rf  = 661;
% ci  = 474;
% cf  = 745;
% f4  = figure(4);
% set(f4, 'position', [112 177 560 420])
% subplot(2,1,1)
% imagesc(imdata0(ri:rf,ci:cf))
% colormap(gray)
% title(['TOP feature in ', fname0])
% axis equal tight off
% xlabel('y (pixels)')
% ylabel('x (pixels)')
% 
% subplot(2,1,2)
% imagesc(imdata(ri:rf,ci:cf))
% colormap(gray)
% title(['TOP feature in ', fname])
% axis equal tight off
% xlabel('y (pixels)')
% ylabel('x (pixels)')
