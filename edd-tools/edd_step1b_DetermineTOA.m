clear all
close all
clc

proot   = '/home/s1b/__eval/projects_parkjs/edd_6bm_2017-1/startup_mar17';

format long

%%%%%%%%%%%%%%%%%
% INPUT
% CHANNEL TO ENERGY CONVERSION RESULT FROM PREVIOUS STEP USING Cd109
%%%%%%%%%%%%%%%%%
%%%
det_id  = 1;
ChToEnergyConversion    =  [0.0347612 0.0375124]; % det1
TOA  = 5.51399;
%%%
% det_id  = 2;
% ChToEnergyConversion    = [0.0341563 -0.00480304]; % det2
% TOA  = 5.24072;
%%%

%%%%%%%%%%%%%%%
% Nominal ExperimentalTOA  = 5.95939; Geometry
TOA0    = 5.5;

%%%%%%%%%%%%%%%
% X-ray emission lines (eV) - XRAY ORANGE BOOK
% Ce Ka1     Ka2     Kb1     La1    La2    Lb1    Lb2    Lg1  Ma1
CeO2_emission_energy        = load('Ce.emission.data');
CeO2_emission_energy(1,:)   = CeO2_emission_energy(1,:)/1000;

%%%%%%%%%%%%%%%
% CeO2 lattice constant
% 5.411651
% fcc
LattParms   = 5.411651;
hkls        = load('fcc.hkls')';
d_hkl       = PlaneSpacings(LattParms, 'cubic', hkls);
lambda_hkl0 = 2.*d_hkl*sind(TOA0/2);
E_hkl0      = Angstrom2keV(lambda_hkl0);
peaks2use   = [2 3 4 5 6 7 9];   %%% USE 4 CeO2 PEAKS TO GET TOA

%%%%%%%%%%%%%%%%%
% USE Ceria diffraction data
pname   = 'CeO2_hv_gv_20170314';
froot   = pname;

fname   = sprintf('%s-%03d-hv.xy', froot, 50);
pfname  = fullfile(proot, pname, fname);

data    = load(pfname);

x   = 1:1:8192;
x   = 1:1:8192;
if det_id == 1
    y   = data(:,1); %%% ChToEnergyConversion    = [0.0347612 0.0375123];
elseif det_id == 2
    y   = data(:,2); %%% ChToEnergyConversion    = [0.0341563 -0.00480304];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_grid  = ChToEnergyConversion(1)*x + ChToEnergyConversion(2);

figure(1)
plot(E_grid, y, 'b.');
hold on
plot(CeO2_emission_energy(1,:), CeO2_emission_energy(2,:), 'g^')
plot(E_hkl0, ones(length(E_hkl0), 1), 'r^')
grid on
xlabel('energy (keV)')
ylabel('counts')

for i = 1:1:length(peaks2use)
    E0  = E_hkl0(peaks2use(i));

    idx1    = find(E_grid < (E0 - 1.5));
    idx2    = find(E_grid < (E0 + 1.5));
    idx1    = idx1(end);
    idx2    = idx2(end);
    
    xdata   = E_grid(idx1:idx2);
    ydata   = y(idx1:idx2);
    
    p0  = [ ...
        max(ydata); ...
        2.5; ...
        0.5; ...
        E0; ...
        0; ...
        0; ...
        ];
    
    pLB = [0; 0; 0; E0-3; -inf; -inf];
    pUB = [inf; inf; 1; E0+3; inf; inf];
    
    p   = lsqcurvefit(@pfunc, p0, xdata, ydata, pLB, pUB);
    
    A(i)    = p(1);
    
    yfit0   = pfunc(p0, xdata);
    yfit    = pfunc(p, xdata);
    
    plot(xdata, ydata, 'm.')
    plot(xdata, yfit0, 'k:')
    plot(xdata, yfit, 'k-')
    
    E_fit(i)    = p(4);
end
lambda_fit  = keV2Angstrom(E_fit);

TOA             = lsqcurvefit(@funcBragg, TOA0, d_hkl(peaks2use), lambda_fit);
lambda_hkl      = 2.*d_hkl*sind(TOA/2);
E_hkl           = Angstrom2keV(lambda_hkl);
pseudo_strain   = E_fit./E_hkl(peaks2use) - 1;

figure(1)
plot(E_hkl, ones(length(E_hkl), 1), 'k^')

figure(2)
plot(2.*d_hkl(peaks2use), lambda_fit, 'ko')
hold on
plot(2.*d_hkl(peaks2use), lambda_hkl(peaks2use), 'g-')
xlabel('2*d_{hkl} (Angstrom)')
ylabel('lambda (Andstrom)')
legend('data points', 'Bragg func fit')

figure(3)
plot(2.*d_hkl(peaks2use), pseudo_strain, 'ko')
xlabel('2*d_{hkl} (Angstrom)')
ylabel('pseudo-strain (-)')

figure(4)
plot(x, log(y), 'b.-');
hold on
plot((CeO2_emission_energy(1,:)-ChToEnergyConversion(2))/ChToEnergyConversion(1), log(CeO2_emission_energy(2,:)), 'g^')
plot((E_hkl-ChToEnergyConversion(2))/ChToEnergyConversion(1), ones(length(E_hkl), 1), 'k^')
grid on
xlabel('channel number')
ylabel('counts')

disp(sprintf('*****************************************************************'))
disp(sprintf('Channel to Energy relathionship'));
disp(sprintf('E (keV) = m * ChannelNumber + b'));
disp(sprintf('m = %f', ChToEnergyConversion(1)));
disp(sprintf('b = %f', ChToEnergyConversion(2)));
disp(sprintf('Channel Number starts from 1'));
disp(sprintf('*****************************************************************'))
disp(sprintf('lambda_hkl = sin (TOA / 2) * (2 * d_hkl)'))
disp(sprintf('Figure 2 should be linear'))
disp(sprintf('TOA = %f degree', TOA));
disp(sprintf('*****************************************************************'))
disp(sprintf('*************** ATTENTION | ACHTUNG | FIGYELEM ******************'))
disp(sprintf('Copy the these commands into following steps.'))
disp(sprintf('TOA  = %g;', TOA))
disp(sprintf('ChToEnergyConversion  = [%g %g];', ChToEnergyConversion(1), ChToEnergyConversion(2)));
