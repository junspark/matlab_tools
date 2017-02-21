clear all
close all
clc

%%% FOR PIXI1
% pname   = 'C:\Users\parkjs\Documents\GitHub\matlab_tools_examples\image-processing-test-images';
% fname   = 'AgBeh_013043.tif';
% 
% pfname  = fullfile(pname, fname);
% 
% raw = imread(pfname);
% % csq = ReadPixirad(pfname, 'nxsq', 900);
% csq = ReadPixirad(pfname);
% 
% figure(1)
% imagesc(log(double(raw)))
% % imagesc(raw)
% axis equal
% 
% figure(2)
% imagesc(log(csq))
% % imagesc(csq)
% axis equal

%%% FOR PIXI2
pname   = '/home/beams/S1IDUSER/mnt/s1c/pixirad_feb17';
% fname   = 'test_001995.tif';
fname   = 'furnace_test_000030.tif';

pfname  = fullfile(pname, fname);

raw = imread(pfname);
csq = ReadPixirad(pfname, 'version', 'pixi2');

figure(1)
imagesc(raw)
% imagesc(raw)
axis equal tight

figure(2)
imagesc(csq)
% imagesc(csq)
axis equal tight

% %%% CORRECTION TABLE FOR PIXI2
fid = fopen('/home/beams/S1IDUSER/mnt/s1a/misc/pixirad2/usb.after_repair/Calibrations/2010.crrm', 'r');
data    = fread(fid, 'float');
fclose(fid);
data    = reshape(data, 1024, 402);

figure(3)
imagesc(data)