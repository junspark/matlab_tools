clear all
% close all
clc

LS_30nm = load('AveLatticeStrain.30nm_norminc.mat');
LS_500nm = load('AveLatticeStrain.500nm_norminc.mat');

figure(10001)
plot((1:1:LS_30nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_30nm.AveLatticeStrain{1}(1,:), 'r^')
hold on
plot((1:1:LS_500nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_500nm.AveLatticeStrain{1}(1,:), 'bv')
title('Cu 111')
xlabel('azimuth (deg)')
ylabel('strain')
grid on
axis([0 350 -3e-3 3e-3])
legend('30 nm' , '500 nm')

figure(10002)
plot((1:1:LS_30nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_30nm.AveLatticeStrain{1}(2,:), 'r^')
hold on
plot((1:1:LS_500nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_500nm.AveLatticeStrain{1}(2,:), 'bv')
title('Cu 200')
xlabel('azimuth (deg)')
ylabel('strain')
grid on
axis([0 350 -3e-3 3e-3])
legend('30 nm' , '500 nm')

figure(10003)
plot((1:1:LS_30nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_30nm.AveLatticeStrain{2}(1,:), 'r^')
hold on
plot((1:1:LS_500nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_500nm.AveLatticeStrain{2}(1,:), 'bv')
title('Nb 110')
xlabel('azimuth (deg)')
ylabel('strain')
grid on
axis([0 350 -3e-3 3e-3])
legend('30 nm' , '500 nm')

figure(10004)
plot((1:1:LS_30nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_30nm.AveLatticeStrain{2}(2,:), 'r^')
hold on
plot((1:1:LS_500nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_500nm.AveLatticeStrain{2}(2,:), 'bv')
title('Nb 200')
xlabel('azimuth (deg)')
ylabel('strain')
grid on
axis([0 350 -3e-3 3e-3])
legend('30 nm' , '500 nm')

figure(10005)
plot((1:1:LS_30nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_30nm.AveLatticeStrain{2}(3,:), 'r^')
hold on
plot((1:1:LS_500nm.XRDIMAGE.CakePrms.bins(1)).*10, LS_500nm.AveLatticeStrain{2}(3,:), 'bv')
title('Nb 211')
xlabel('azimuth (deg)')
ylabel('strain')
grid on
axis([0 350 -3e-3 3e-3])
legend('30 nm' , '500 nm')