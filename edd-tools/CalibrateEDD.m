%%%
% TEMPLATE FOR EDD DETECTEOR CALIBRAION
%%%
clear all
close all
clc

pname = 'C:\Users\parkjs\Documents\GitHub\matlab_tools\edd-tools\calibration-examples';
fname ='cal_7d_pb_100s';
pfname = fullfile(pname, fname);

[x,y] = ReadEDDData(pfname);

figure(1)
plot(x,y)
xlabel('MCA Channel Number')
ylabel('counts (arb. units)')
grid on
hold on

Bounds  = [ ...
    328 379; ...
    379 400; ...
    400 423; ...
    472 526; ...
    526 566; ...
    566 581; ...
    581 605; ...
    621 655; ...
    865 925; ...
    1015 1075; ...
    1084 1120; ...
    1248 1294; ...
    1360 1400; ...
    1400 1445; ...
    1525 1585; ...
    1620 1675; ...
    ];

xfit    = zeros(size(Bounds,1),1);
for i = 1:1:size(Bounds,1)
    xdata   = x(Bounds(i,1):Bounds(i,2));
    ydata   = y(Bounds(i,1):Bounds(i,2));
    
    [A0, x0]    = max(ydata);
    x0  = x0 + Bounds(i,1);
    
    p0  = [ ...
        A0; ...
        3.5; ...
        0.5; ...
        x0; ...
        0; ...
        0; ...
        ];    
    
    p	= lsqcurvefit(@pfunc, p0, xdata, ydata);
    xfit(i,1)   = p(4);
    
    yfit0   = pfunc(p0, xdata);
    yfit    = pfunc(p, xdata);
    
    figure(1)
    plot(xdata, ydata, 'r.')
    plot(xdata, yfit0, 'g.-')
    plot(xdata, yfit, 'k.-')
end

th      = 3;
a0      = 5.411102;
hkls    = load('fcc.hkls');
hkls2   = sqrt(sum(hkls.*hkls,2));
d0_hkls = a0./hkls2;
lambda  = (2*sind(th)).*d0_hkls;
Energy  = Angstrom2keV(lambda);

figure(2)
plot(hkls2, Energy, 'k.')
xlabel('sqrt(hkls)')
ylabel('energy (keV)')
grid on

xCeO2   = [xfit(5); xfit(8:end)];
Energy  = Energy(1:10);

figure(3)
plot(hkls2(1:10), xCeO2./hkls2(1:10), 'b^:')
hold on
plot(hkls2(1:10), Energy./hkls2(1:10), 'k.-') 

% xCeO2   = xfit;
% Energy  = Energy;
% 
% figure(3)
% plot(hkls2(1:16), xCeO2./hkls2(1:16), 'b^:')
% hold on
% plot(hkls2, Energy./hkls2, 'k.-') 

pCalibration    = polyfit(xCeO2, Energy, 1)
EnergyFit       = polyval(pCalibration, xCeO2);

figure(4)
plot(xCeO2, Energy, 'b.')
hold on
plot(xCeO2, EnergyFit, 'k-') 