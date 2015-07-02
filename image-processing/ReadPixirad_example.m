clear all
close all
clc

pname   = 'C:\Users\parkjs\Desktop\temp';
fname   = 'pixiradtest_248926.tiff';

pfname  = fullfile(pname, fname);

raw = imread(pfname);
csq = ReadPixirad(pfname, 'nxsq', 900);

figure(1)
imagesc(log(double(raw)))
% imagesc(raw)
axis equal

figure(2)
imagesc(log(csq))
% imagesc(csq)
axis equal
return




