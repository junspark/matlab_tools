clear all
close all
clc

pname   = './';
fname   = 'pixiradtest_248926.tiff';

pfname  = fullfile(pname, fname);

raw = imread(pfname);
% csq = ReadPixirad(pfname, 'nxsq', 900);
csq = ReadPixirad(pfname);

figure(1)
imagesc(log(double(raw)))
% imagesc(raw)
axis equal

figure(2)
imagesc(log(csq))
% imagesc(csq)
axis equal
return




