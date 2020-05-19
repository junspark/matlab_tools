clear all
close all
clc

addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

pname_img               = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\data\dexela';
pname_img_corrected0    = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\data\dexela\corrected';

% froot   = 'am316_ss_id105_1';
% imnum   = 324:443;

% froot   = 'am316_ss_id105_2';
% imnum   = 444:563;

% froot   = 'am316_ss_id105_ht_1050c';
% imnum   = 564:923;

% froot   = 'am316_ss_id206_2';
% imnum   = 204:323;

froot   = 'am316_ss_id206';
imnum   = 84:203;

pname_img_corrected = fullfile(pname_img_corrected0, froot);
if ~isdir(pname_img_corrected)
    mkdir(pname_img_corrected)
end
for i = 1:1:length(imnum)
    fname_img   = sprintf('%s_%06d.tif', froot, imnum(i));
    pfname_img  = fullfile(pname_img, fname_img);
    
    imdata  = imread(pfname_img);
    
    [nnn, mmm]  = size(imdata);
    bkg_synthetic   = median(imdata, 1);
    bkg_synthetic   = repmat(bkg_synthetic, nnn, 1);
    
    imcorrected = flipud(imdata - bkg_synthetic);
    
    fname_img_out   = sprintf('%s_%06d.bkg_corrected.tif', froot, imnum(i));
    pfname_img_out  = fullfile(pname_img_corrected, fname_img_out);
    WriteDexela(pfname_img_out, imcorrected);
end

figure,
imagesc(imdata)
axis equal tight
caxis([-50 1000])

figure,
imagesc(imcorrected)
axis equal tight
caxis([-50 1000])

% imageinfo(pfname_img)
% imageinfo(fname_img_out)