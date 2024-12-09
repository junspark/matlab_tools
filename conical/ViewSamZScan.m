clear all
close all
clc

%%% analyse samZ scans
% load background
pname   = '~/mnt/s1c/brown_aug17/ge/Sample_D3';
froot   = 'dark2s_';
fnumber = 2062;

fname   = sprintf('%s%06d.ge3', froot, fnumber);
pfname  = fullfile(pname, fname);
bg  = NreadGE(pfname, 1);

froot = 'cs_testing200_';
samZ = linspace(-6, 6, 61);
imno = 2001:2061;

for mm = length(imno):-1:1
    fnumber = imno(mm);
    fname   = sprintf('%s%06d.ge3', froot, fnumber);
    pfname  = fullfile(pname, fname);
    
    im = NreadGE(pfname, 1);
    img = rot90(im-bg);
    mm_sum(mm)  = sum(img(:));
    
%     figure(1)
%     imagesc(img,[-10 500]); colorbar
%     axis square
%     title(imno(mm))
end

plot(samZ, mm_sum)