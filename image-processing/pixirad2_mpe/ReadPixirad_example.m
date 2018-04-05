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
pfname  = '/home/beams/S1IDUSER/mnt/s1c/birkedal_feb17/pixirad/glassyC_1s_000181.tif';
csq = ReadPixirad(pfname, 'version', 'pixi2');

figure(1)
imagesc(log(csq))
