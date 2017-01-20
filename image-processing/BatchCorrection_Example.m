clear all
close all
clc

path_bkg    = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dm_batchcorr';
bkg_num     = 1;
root_bkg    = 'dark_';
genum       = 2;

path_image  = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dm_batchcorr';
root_image  = 'Ni_test_';

path_output = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dm_batchcorr';

BatchCorrection(path_bkg, bkg_num, root_bkg, ...
    path_image, root_image, ...
    genum, path_output, ...
    'CorrectAllImages', false, ...
    'lo', 10000, 'hi', 10001, ...
    'OutAllFrames', false, ...
    'DisplayFrames', false)