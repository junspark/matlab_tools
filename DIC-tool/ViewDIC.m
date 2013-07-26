clear all
close all
clc

pname   = 'V:\Bieler_July13\DIC';
fname0  = 'DIC_00045.tif'; % Initial state
pfname0 = fullfile(pname, fname0);
imdata0 = imread(pfname0);

fname   = 'DIC_00506.tif';      %%%%%%%%%%%%%%%%
pfname  = fullfile(pname, fname);
imdata  = imread(pfname);

% imdata0 = imdata0./mean(imdata0(:));
% imdata  = imdata./mean(imdata(:));

% imdata0 = double(imdata0)*hamming(length(imdata0));
% imdata  = double(imdata)*hamming(length(imdata));

filter0 = blackman(size(imdata0,1),'symmetric'); % setup window in one direction - I chose Blackman window - could use Hamming or others
filter  = blackman(size(imdata,2),'symmetric'); % setup window in second direction
filtimg = (filter0*filter'); % make 2D windowing function that is sized to match image
imdata0 = double(imdata0).*filtimg; % apply windowing function to image
imdata  = double(imdata).*filtimg;

FFT0    = fft2(imdata0); % 2d FFT
FFT     = conj(fft2(imdata));
FFTR    = FFT0.*FFT;
magFFTR = abs(FFTR);
FFTRN   = (FFTR./magFFTR);

result = ifft2(double(FFTRN));

figure;
colormap('gray'); % choose colormap for plotting
imagesc(result);

figure;
colormap(jet);
mesh(result); 
return
ri  = 593;
rf  = 900;
ci  = 20;
cf  = 870;
x   = (ci:1:cf) - ci;
x   = x*0.0019;

profile0    = imdata0(ri:rf,ci:cf);
profile0    = sum(profile0,1);

profile     = imdata(ri:rf,ci:cf);
profile     = sum(profile,1);

f1  = figure(1);
% set(f1, 'position', [680 678 560 420])
subplot(2,2,1)
imagesc(imdata0)
title(fname0)
axis equal tight off
xlabel('y')
ylabel('x')

subplot(2,2,2)
imagesc(imdata)
title(fname)
axis equal tight off
xlabel('y')
ylabel('x')

subplot(2,2,3)
imagesc(imdata0 - imdata)
title('difference image')
axis equal tight off
xlabel('y')
ylabel('x')

subplot(2,2,4)
plot(x, profile0, 'rs-')
hold on
plot(x, profile, 'bo-')
xlabel('y (mm)')
ylabel('edge (mm)')

ri  = 324;
rf  = 600;
ci  = 2152;
cf  = 2440;
f2  = figure(2);
set(f2, 'position', [1249 178 560 420])
subplot(2,1,1)
imagesc(imdata0(ri:rf,ci:cf))
title(['BOTTOM ROI - feature in ', fname0])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

subplot(2,1,2)
imagesc(imdata(ri:rf,ci:cf))
title(['BOTTOM ROI - feature in ', fname])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

ri  = 675;
rf  = 723;
ci  = 1369;
cf  = 1555;
f3  = figure(3);
set(f3, 'position', [681 178 560 420])
subplot(2,1,1)
imagesc(imdata0(ri:rf,ci:cf))
title(['BEAM ROI - feature in ', fname0])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

subplot(2,1,2)
imagesc(imdata(ri:rf,ci:cf))
title(['BEAM ROI - feature in ', fname])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

ri  = 541;
rf  = 661;
ci  = 474;
cf  = 745;
f4  = figure(4);
set(f4, 'position', [112 177 560 420])
subplot(2,1,1)
imagesc(imdata0(ri:rf,ci:cf))
title(['TOP feature in ', fname0])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')

subplot(2,1,2)
imagesc(imdata(ri:rf,ci:cf))
title(['TOP feature in ', fname])
axis equal tight off
xlabel('y (pixels)')
ylabel('x (pixels)')
