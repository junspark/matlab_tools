clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCC Fe
latticeParms_bcc    = 2.82 ;                        % IN Angstrom
hkls_bcc            = load('bcc.hkls');
d_hkl_bcc           = PlaneSpacings(latticeParms_bcc, 'cubic', hkls_bcc');

% FCC Fe
latticeParms_fcc    = 3.515;                        % IN Angstrom
hkls_fcc            = load('fcc.hkls');
d_hkl_fcc           = PlaneSpacings(latticeParms_fcc, 'cubic', hkls_fcc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HORZ
ChToEnergyConversion = [0.0926   -0.0754];
TOA  = 6.9678;

% % %%% VERT
% % ChToEnergyConversion = [0.0976   -0.2267];
% % TOA =  6.6375;

lambda_hkl_bcc  = 2.*d_hkl_bcc*sind(TOA/2);
lambda_hkl_fcc  = 2.*d_hkl_fcc*sind(TOA/2);

E_hkl_bcc   = Angstrom2keV(lambda_hkl_bcc);
E_hkl_fcc   = Angstrom2keV(lambda_hkl_fcc);

pname_data  = 'C:\Users\parkjs\Box Sync\Projects\Operations\UserSupport\mach_201502_edd\Lap7';
dirlist = dir([pname_data, '/horz*']);
data    = zeros(2048,1);
dataStack   = [];
for i = 1:1:length(dirlist)
    pfname  = fullfile(pname_data, dirlist(i).name);
    datai   = load(pfname);
    data    = data + datai;
    dataStack   = [dataStack datai];
end

x   = 1:1:2048;
E   = ChToEnergyConversion(1)*x + ChToEnergyConversion(2);

figure,
plot(E, data, 'k.')
hold on
plot(E_hkl_bcc, ones(length(E_hkl_bcc), 1), 'b^')
% plot(E_hkl_fcc, ones(length(E_hkl_fcc), 1), 'r^')
grid on
xlabel('Energy (keV)')
ylabel('counts')

figure,
imagesc(log(dataStack))
set(gca, 'YTickLabel', num2str(round(str2num(get(gca, 'YTickLabel'))*ChToEnergyConversion(1) + ChToEnergyConversion(2))))
title('horz')
xlabel('data point number')
ylabel('Channel Number => Energy')