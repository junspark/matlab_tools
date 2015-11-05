clear all
close all
clc
load wscub4x
frmesh  = wscub4x.frmesh;

% Load MIDAS results
pfname  = 'Grains.csv';
Grains  = parseGrainData(pfname, frmesh.symmetries);
numpts  = length(Grains);
wts     = ones(numpts, 1);
quat    = [Grains.quat];

odf = CalcOdfFromAggregate(frmesh, quat, wts, 'PlotOdf', 'on', 'std', 10);
