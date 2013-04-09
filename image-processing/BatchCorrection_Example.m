clear all
close all
clc

path_bkg    = 'C:\Users\s1iduser\Documents\parkjs\batchcorr_example';
bkg_num     = 1;
root_bkg    = 'dark_';
genum       = 2;

path_image  = 'C:\Users\s1iduser\Documents\parkjs\batchcorr_example';
root_image  = 'Narwhale_3_';

path_output = 'C:\Users\s1iduser\Documents\parkjs\batchcorr_example\corrected';

BatchCorrection(path_bkg, bkg_num, root_bkg, ...
    path_image, root_image, ...
    genum, path_output, ...
    'CorrectAllImages', false, ...
    'lo', 1966, 'hi', 1966, ...
    'OutAllFrames', true, ...
    'DisplayFrames', false)