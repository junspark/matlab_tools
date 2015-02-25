clear all
close all
clc

format long
%%%%%%%%%%%%%%%
% Cd109 radioactive decay data (keV) - lbl.gov
% E  88.045 2.634 2.806 2.978 2.984 3.15 3.203 3.234 3.256 3.348 3.520 3.743 3.750 21.708 21.990 22.163 24.912 24.943 25.144 25.455 25.511
% I    3.61 0.183 0.508 4.57 2.64 0.144 0.226 0.0305 0.589 0.284 0.0277 0.045 0.00122 29.5 55.7 4.79 9.23 0.0673 2.31 0.487
Cd109_energy	= load('Cd109.decay.data');

%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%
% USE Cd109 decay data
pname_Cd109_spec    = '.\calibration-examples\mach_feb15_calibration\horizontal';
fname_Cd109_spec    = 'Cd109_gain1p70_maxE40real';
pfname_Cd109_spec   = fullfile(pname_Cd109_spec, fname_Cd109_spec);
[x, y]  = ReadEDDData(pfname_Cd109_spec, 'IDLFile', 1);

%%%%%%%%%%%%%%%%%
% INITIAL GUESS FOR DECAY PEAKS
x0  = [240; 271; 952];
y0  = [13570; 3079; 591].*20;
g0  = ones(length(x0), 1).*20;
n0  = ones(length(x0), 1).*0.5;
p0  = [y0 g0 n0 x0]';
p0  = p0(:);
p0  = [p0; 0; 1];
pLB = [0; 0; 0; x0(1)-10; 0; 0; 0; x0(1)-10; 0; 0; 0; x0(1)-10; -inf; -inf];
pUB = [inf; inf; 1; x0(1)+10; inf; inf; 1; x0(2)+10; inf; inf; 1; x0(3)+10; inf; inf];
%%%%%%%%%%%%%%%%%

p       = lsqcurvefit(@pfunc, p0, x, y, pLB, pUB);

yfit0   = pfunc(p0, x);
yfit    = pfunc(p, x);

figure(1)
plot(x, y, 'r.');
hold on
plot(x, yfit0, 'g:')
plot(x, yfit, 'k-');
grid on
xlabel('channel number')
ylabel('counts')    
set(g1, 'Position', [0.13 0.345238 0.775 0.579762])

g2  = subplot(2,1,2);
plot(x, (y-yfit0)./y, 'g:');
hold on
plot(x, (y-yfit)./y, 'k-');
grid on
xlabel('channel number')
ylabel('difference ratio')    
set(g2, 'Position', [0.130 0.115 0.775 0.118])

ChToEnergyEquivalence   = [ ...
    p(4) sum(Cd109_energy(14:16,1).*Cd109_energy(14:16,2))./sum(Cd109_energy(14:16,2)); ...
    p(8) sum(Cd109_energy(17:21,1).*Cd109_energy(17:21,2))./sum(Cd109_energy(17:21,2)); ...
    p(12) Cd109_energy(1); ...
    ];
ChToEnergyConversion    = polyfit(ChToEnergyEquivalence(:,1), ChToEnergyEquivalence(:,2) , 1)
