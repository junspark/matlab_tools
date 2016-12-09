clear all
close all
clc

format long
%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%
% USE Cd109 decay data
pname_Cd109_spec    = '/home/s1b/__eval/edd_mach_dec16/cal_Co57_20161207/';
fname_Cd109_spec    = 'cal_Co57_20161207-001-hv.xy';

pfname_Cd109_spec   = fullfile(pname_Cd109_spec, fname_Cd109_spec);
data    = load(pfname_Cd109_spec);
x   = 1:1:8192;
yv  = data(:,1);
yh  = data(:,2);

y   = yv; % ChToEnergyConversion    = [0.0347817 0.0201449]; v
y   = yh; % ChToEnergyConversion    = [0.034798 0.0112459];; h

%%%%%%%%%%%%%%%
% Cd109 radioactive decay data (keV) - lbl.gov
% E  88.045 2.634 2.806 2.978 2.984 3.15 3.203 3.234 3.256 3.348 3.520 3.743 3.750 21.708 21.990 22.163 24.912 24.943 25.144 25.455 25.511
% I    3.61 0.183 0.508 4.57 2.64 0.144 0.226 0.0305 0.589 0.284 0.0277 0.045 0.00122 29.5 55.7 4.79 9.23 0.0673 2.31 0.487
Co57_energy     = load('Co57.decay.data');

%%%%%%%%%%%%%%%%%
% INITIAL GUESS FOR DECAY PEAKS
x0  = [413; 3509; 3924;];  % v
y0  = [5000; 90570; 7000;].*10;
g0  = ones(length(x0), 1).*10;
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

Rwp     = ErrorRwp(y, yfit);
Re      = ErrorRe(y, yfit);
Rp      = ErrorRp(y, yfit);

figure(1)
axes('Position', [0.130 0.345238 0.775 0.579762]);
hold on
plot(x, y, 'r.');
plot(x, yfit0, 'g:')
plot(x, yfit, 'k-');
grid on
xlabel('channel number')
ylabel('counts')    
axis([x(1) x(end) min(y) max(y)])

axes('Position', [0.130 0.115 0.775 0.118]);
hold on
plot(x, (y-yfit0)./y, 'g:');
plot(x, (y-yfit)./y, 'k-');
grid on
xlabel('channel number')
ylabel('difference ratio')    
axis([x(1) x(end) min((y-yfit0)./y) max((y-yfit0)./y)])

ChToEnergyEquivalence   = [ ...
    p(4:4:12)  Co57_energy(:,1); ...
    ];

ChToEnergyConversion    = polyfit(ChToEnergyEquivalence(:,1), ChToEnergyEquivalence(:,2) , 1);

disp(sprintf('Channel to Energy relathionship'));
disp(sprintf('E (keV) = m * ChannelNumber + b'));
disp(sprintf('m = %f', ChToEnergyConversion(1)));
disp(sprintf('b = %f', ChToEnergyConversion(2)));
disp(sprintf('Channel Number starts from 1'));
disp(sprintf('************* ATTENTION | ACHTUNG | FIGYELEM *****************'))
disp(sprintf('Copy the following command into LINE 8 on step1_DetermineTOA.m'))
disp(sprintf('ChToEnergyConversion    = [%g %g];', ChToEnergyConversion(1), ChToEnergyConversion(2)));
