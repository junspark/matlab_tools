clear all
close all
clc
%%%%%%%%%%%%%%%%%
% INPUT
% CHANNEL TO ENERGY CONVERSION RESULT FROM PREVIOUS STEP USING Cd109
%%%%%%%%%%%%%%%%%
ChToEnergyConversion    = [0.03928 0.0754175];

%%%%%%%%%%%%%%%
% Nominal Experimental Geometry
TOA0    = 5.468202;

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
peaks2use   = [2 3 4];   %%% USE 4 CeO2 PEAKS TO GET TOA

%%%%%%%%%%%%%%%%%
% USE Ceria diffraction data
ngrid   = 31;
fstem   = 'cal1_ceria_gv';
pname	= './cal1_ceria_gv';
for i = 1:1:ngrid
    fname_CeO2_spec    = sprintf('%s-%03d.xy', fstem, i);
    pfname_CeO2_spec   = fullfile(pname, fname_CeO2_spec);
    [x, y]  = ReadEDDData(pfname_CeO2_spec, 'SpecFile', 1);
end