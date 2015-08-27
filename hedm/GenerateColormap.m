% NEED TO BE IN 
% /net/s1dserv/export/s1b/__eval/mtex-4.0.23
% run "startup_mtex" to run these commands to generate colormaps

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v
% GENERATE IPF COLORMAP USING MTEX
% CHECK IN ebsdColorbar.m
cs  = crystalSymmetry('m-3m');
ss  = specimenSymmetry('-1');
cc  = get_option('antipodal','colorcoding','ipdfHKL');

[minTheta,maxTheta,minRho,maxRho,v] = getFundamentalRegionPF(cs, 'antipodal');
h   = S2Grid('PLOT', 'minTheta', minTheta, 'maxTheta', maxTheta,...
    'minRho', minRho, 'maxRho', maxRho, 'RESTRICT2MINMAX', 'resolution', 1*degree, 'antipodal');
v   = vector3d(h);
x   = getx(v); x = x(:);
y   = gety(v); y = y(:);
z   = getz(v); z = z(:);
d   = orientation2color(h,cc,cs,'antipodal');
save('coloring_scheme_cubic.mat', 'x', 'y', 'z', 'd')
return