clear all
close all
clc

%%%%%%%%%%%%%%%%%
% INPUT
% CHANNEL TO ENERGY CONVERSION RESULT FROM PREVIOUS STEP USING Cd109
%%%%%%%%%%%%%%%%%
TOA  = 6.87555;
ChToEnergyConversion  = [0.096902 -0.107915];

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
lambda_hkl	= 2.*d_hkl*sind(TOA/2);
E_hkl       = Angstrom2keV(lambda_hkl);
peaks2use   = [3 4 5 6 7 8 9 10];   %%% USE 4 CeO2 PEAKS TO GET TOA

%%%%%%%%%%%%%%%%%
% USE Ceria diffraction data
pname_CeO2_spec    = '.\calibration-examples\edd_jun15_calibration';    % V
fname_CeO2_spec    = 'ceria_6june2015_V_maxe40_300s';             % V
pfname_CeO2_spec   = fullfile(pname_CeO2_spec, fname_CeO2_spec);
[x, y]  = ReadEDDData(pfname_CeO2_spec, 'IDLFile', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_grid  = ChToEnergyConversion(1)*x + ChToEnergyConversion(2);

figure(1)
plot(E_grid, y, 'b.');
hold on
plot(CeO2_emission_energy(1,:), CeO2_emission_energy(2,:), 'g^')
plot(E_hkl, ones(length(E_hkl), 1), 'r^')
grid on
xlabel('energy (keV)')
ylabel('counts')

pkfitting_pars.pfunc_type   = 'splitpseudovoigt';
pkfitting_pars.pbkg_order   = 2;
for i = 1:1:length(peaks2use)
    E   = E_hkl(peaks2use(i));

    idx1    = find(E_grid < (E - 1.5));
    idx2    = find(E_grid < (E + 1.5));
    idx1    = idx1(end);
    idx2    = idx2(end);
    
    xdata   = E_grid(idx1:idx2);
    ydata   = y(idx1:idx2);
    
    p0  = [ ...
        max(ydata); ...
        0.5; ...
        0.5; ...
        0.01; ...
        0.01; ...
        E; ...
        0; ...
        0; ...
        ];
    pLB = [ ...
        0; ...
        0; ...
        0; ...
        0; ...
        0; ...
        E - 1.5; ...
        -inf; ...
        -inf; ...
        ];
    pUB = [ ...
        inf; ...
        inf; ...
        inf; ...
        1; ...
        1; ...
        E + 1.5; ...
        inf; ...
        inf; ...
        ];
    
    pkfitting_pars.xdata    = xdata;
    [p, rn(i), ~, ef(i)]    = lsqcurvefit(@pfunc_switch, p0, pkfitting_pars, ydata, pLB, pUB);
    
    yfit0   = pfunc_switch(p0, pkfitting_pars);
    yfit    = pfunc_switch(p, pkfitting_pars);
    
    plot(xdata, ydata, 'm.')
    plot(xdata, yfit0, 'k:')
    plot(xdata, yfit, 'k-')
    
    E_fit(i)        = p(6);
    lambda_fit(i)   = keV2Angstrom(E_fit(i));
end
pseudo_strain   = E_fit./E_hkl(peaks2use) - 1

figure(1)
plot(E_fit, ones(length(E_fit), 1), 'k^')

figure(2)
subplot(2,1,1)
plot(2.*d_hkl(peaks2use), lambda_fit, 'ko')
hold on
plot(2.*d_hkl(peaks2use), lambda_hkl(peaks2use), 'g-')
xlabel('2*d_{hkl} (Angstrom)')
ylabel('lambda (Andstrom)')
legend('data points', 'Bragg func fit')
grid on

subplot(2,1,2)
plot(2.*d_hkl(peaks2use), pseudo_strain, 'ko-')
grid on
xlabel('2*d_{hkl} (Angstrom)')
ylabel('pseudo-strain (-)')

figure(3)
plot((E_hkl(peaks2use) - ChToEnergyConversion(2))/ChToEnergyConversion(1), pseudo_strain, 'ko-')
grid on
xlabel('MCA channel number (-)')
ylabel('pseudo-strain (-)')