clear all
close all
clc

%%% Pb
a0  = 4.9508;
hkls    = load('fcc.hkls');
hkls2   = sqrt(sum(hkls.*hkls,2));
d0_hkls = a0./hkls2;

pname   = 'W:\__eval\park_jul14_edd\data\calibration_filed_for_horizontal';
fname   ='cal_7d_pb_100s';
pfname  = fullfile(pname, fname);
ydata   = dlmread(pfname, '', 12, 0);
xdata   = 1:1:length(ydata);

pEDD = [
    0.0498
   -0.0330
    3.5000];
xdata_d_hkl = EDDModel(pEDD,xdata);

figure,
plot(xdata, ydata)

figure,
plot(pEDD(1)*xdata + pEDD(2), ydata)
return
figure,
plot(xdata_d_hkl(200:end), ydata(200:end))
hold on
plot(d0_hkls, 10, 'r^')

%%% Pb
a0  = 4.9508;
hkls    = load('fcc.hkls');
hkls2   = sqrt(sum(hkls.*hkls,2));
d0_hkls = a0./hkls2;

pname   = 'W:\__eval\park_jul14_edd\data\calibration_for_vertical_plane';
fname   ='calPbV_july12.dat';
pfname  = fullfile(pname, fname);
ydata   = dlmread(pfname, '', 97, 0);
xdata   = 1:1:length(ydata);

pEDD = [
    0.0498
   -0.0330
    3.5000];
xdata_d_hkl = EDDModel(pEDD,xdata);


figure,
plot(xdata, ydata)

% figure,
% plot(pEDD(1)*xdata + pEDD(2), ydata)
% 
% figure,
% plot(xdata_d_hkl(200:end), ydata(200:end))
% hold on
% plot(d0_hkls, 10, 'r^')