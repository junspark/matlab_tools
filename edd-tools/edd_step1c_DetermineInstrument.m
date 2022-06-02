clear all
close all
clc

addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% EXAMPLE DATA SET1
% pfname_metadata = './AFRL_Elmer_oct21/exposure_record_full.par';
% proot           = './AFRL_Elmer_oct21';
% froot           = 'ceria_03nov21_gv3_100x100_preci_0';
% 
% % %%%
% % det_id  = 1;
% % ChToEnergyConversion    =  [0.026066103116357 0.205383496511323]; % det1
% % TOA  = 6.443987906709726;
% % %%%
% 
% %%%
% det_id  = 2;
% ChToEnergyConversion    = [0.026531689524013 0.023759670042594]; % det2
% TOA  = 6.942735993911892;
% %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% EXAMPLE DATA SET2
pfname_metadata = './thwang_feb22/exposure_record_full.par';
proot           = './thwang_feb22';
froot           = 'cali_ceria_200x200_gv_0210';

% %%%
% det_id  = 1;
% ChToEnergyConversion    =  [0.026078802888531 0.190064516497945]; % det1
% TOA  = 6.361913075256187;
% %%%

%%%
det_id  = 2;
ChToEnergyConversion    = [0.026501417803144 0.036784574090575]; % det2
TOA  = 6.575727211674684;
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metadata    = ReadSpecParFile(pfname_metadata, 'Version', 'mpe_standard');
pfroot      = fullfile(proot, froot, sprintf('%s.xy', froot));
edd         = readedd5_6bm(pfroot);

format long

%%%%%%%%%%%%%%%
% Nominal Experimental Geometry
TOA0    = TOA;

%%%%%%%%%%%%%%%
% X-ray emission lines (eV) - XRAY ORANGE BOOK
% Ce Ka1     Ka2     Kb1     La1    La2    Lb1    Lb2    Lg1  Ma1
CeO2_emission_energy        = load('Ce.emission.data');
CeO2_emission_energy(1,:)   = CeO2_emission_energy(1,:)/1000;

%%%%%%%%%%%%%%%
% CeO2 lattice constant
% fcc
LattParms   = 5.41165;   %%% 2007 SRM
% LattParms   = 5.41153;   %%% 2017 SRM
hkls        = load('fcc.hkls')';
d_hkl       = PlaneSpacings(LattParms, 'cubic', hkls);
lambda_hkl0 = 2.*d_hkl*sind(TOA0/2);
E_hkl0      = Angstrom2keV(lambda_hkl0);

% peaks2use   = [3 4 5 6 7 8 9];
peaks2use   = [3 4 6 8 9];

%%%%%%%%%%%%%%%%%
% DETERMINE GV COM AND IDENTIFY CeO2 DATA TO LOOK (CAN SUM THE ENTIRE GV
% SCAN AND WORK WITH THAT TOO)
use_center_pattern  = true;

x   = 1:1:8192;
y   = edd.data{det_id};

if use_center_pattern
    ysum       = sum(y,2);
    pts_grid    = edd.motorpos';
    
    p0  = [ max(ysum); std(pts_grid);   0.5; mean(pts_grid); 0; mean(ysum); ];
    pLB = [          0;   0;    0;    min(pts_grid); -inf; -inf; ];
    pUB = [        inf; inf;    1;    max(pts_grid);  inf;  inf; ];
    
    pkfitting_pars.pfunc_type   = 'pseudovoigt';
    pkfitting_pars.pbkg_order   = 2;
    pkfitting_pars.xdata        = pts_grid;
    
    [p, rnv, ~, efv]    = lsqcurvefit(@pfunc_switch, ...
        p0, pkfitting_pars, ysum, pLB, pUB);
    yfit0   = pfunc_switch(p0, pkfitting_pars);
    yfit    = pfunc_switch(p, pkfitting_pars);
    disp(sprintf('GV center   : %f', p(4)));
    disp(sprintf('GV fwhm     : %f', p(2)));
    
    figure(1);
    plot(pts_grid, ysum, 'ko')
    hold on
    plot(pts_grid, yfit0, 'r--')
    plot(pts_grid, yfit, 'b-')
    xlabel(sprintf('%s [mm]', edd.motorname))
    ylabel('summed intensity [cts]')
    
    idx1    = find(pts_grid > p(4));
    idx2    = find(pts_grid < p(4));
    
    idx1    = idx1(1);
    idx2    = idx2(end);
    idx     = [idx1 idx2];
    
    y_for_cal	= (p(4) - pts_grid(idx(1)))/(pts_grid(idx(2)) - pts_grid(idx(1))).*(y(idx(2),:) - y(idx(1),:)) + y(idx(1),:);
else
    y_for_cal	= mean(y,1)';
end

figure(2)
plot(x, y_for_cal, 'k.')
if use_center_pattern
    hold on
    plot(x, y(idx(1),:), 'r.')
    plot(x, y(idx(2),:), 'b.')
end
xlabel('Channel number [-]')
ylabel('Intensity [cts]')
title('data to determine detector parameters')
grid on

%%% CONVERT CHANNEL NUMBER TO E
E_grid  = polyval(ChToEnergyConversion, x);

figure(3)
plot(E_grid, y_for_cal, 'k.');
hold on
plot(CeO2_emission_energy(1,:), CeO2_emission_energy(2,:), 'g^')
plot(E_hkl0, ones(length(E_hkl0), 1), 'r^')
grid on
xlabel('energy [keV]')
ylabel('counts [-]')

for i = 1:1:length(peaks2use)
    E0  = E_hkl0(peaks2use(i));

    idx1    = find(E_grid < (E0 - 1.5));
    idx2    = find(E_grid < (E0 + 1.5));
    idx1    = idx1(end);
    idx2    = idx2(end);
    
    xdata   = E_grid(idx1:idx2)';
    ydata   = y_for_cal(idx1:idx2)';
    
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
    
    p   = lsqcurvefit(@pfunc, p0, xdata(:), ydata(:), pLB, pUB);
    
    A(i)    = p(1);
    
    yfit0   = pfunc(p0, xdata);
    yfit    = pfunc(p, xdata);
    
    figure(3)
    plot(xdata, ydata, 'm.')
    plot(xdata, yfit0, 'k:')
    plot(xdata, yfit, 'k-')
    
    E_fit(i)    = p(4);
end
lambda_fit      = keV2Angstrom(E_fit);

TOA             = lsqcurvefit(@funcBragg, TOA0, d_hkl(peaks2use), lambda_fit);
lambda_hkl      = 2.*d_hkl*sind(TOA/2);
E_hkl           = Angstrom2keV(lambda_hkl);
pseudo_strain   = E_fit./E_hkl(peaks2use) - 1;

figure(3)
plot(E_hkl, ones(length(E_hkl), 1), 'k^')

figure(4)
plot(2.*d_hkl(peaks2use), lambda_fit, 'ko')
hold on
plot(2.*d_hkl(peaks2use), lambda_hkl(peaks2use), 'g-')
xlabel('2*d_{hkl} (Angstrom)')
ylabel('lambda (Andstrom)')
legend('data points', 'Bragg func fit')
grid on

% figure(3)
% plot(2.*d_hkl(peaks2use), pseudo_strain, 'ko')
% xlabel('2*d_{hkl} (Angstrom)')
% ylabel('pseudo-strain (-)')

% figure(4)
% plot(x, log(y), 'b.-');
% hold on
% plot((CeO2_emission_energy(1,:)-ChToEnergyConversion(2))/ChToEnergyConversion(1), log(CeO2_emission_energy(2,:)), 'g^')
% plot((E_hkl-ChToEnergyConversion(2))/ChToEnergyConversion(1), ones(length(E_hkl), 1), 'k^')
% grid on
% xlabel('channel number')
% ylabel('counts')

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
disp(sprintf('******************** TOA OPTIMIZATION ONLY **********************'))
disp(sprintf('Copy the these commands into following steps.'))
disp(sprintf('TOA  = %g;', TOA))
disp(sprintf('ChToEnergyConversion  = [%g %g];', ChToEnergyConversion(1), ChToEnergyConversion(2)));
% return

%%% DATA TO PASS TO THE LSQNONLIN OPERATION
data.x  = x(:);
data.y  = y_for_cal(:);
data.peaks2use  = peaks2use;
data.d_hkl      = d_hkl;
save('obj_data.mat', 'data')

p0      = [TOA ChToEnergyConversion];
pLB     = [TOA-0.2 ChToEnergyConversion-inf];
pUB     = [TOA+0.2 ChToEnergyConversion+inf];
pseudo_strain0  = funcEDDInstrObjective(p0);
[p,RESNORM,pseudo_strain,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(@funcEDDInstrObjective, ...
    p0, pLB, pUB);

disp(sprintf('*****************************************************************'))
disp(sprintf('**************** AFTER INSTRUMENT CALIBRATION *******************'))
disp(sprintf('Copy the these commands into following steps.'))
disp(sprintf('TOA  = %g;', p(1)))
disp(sprintf('ChToEnergyConversion  = [%g %g];', p(2), p(3)));

figure,
plot(d_hkl(peaks2use), pseudo_strain0, 'bo-')
hold on
plot(d_hkl(peaks2use), pseudo_strain, 'ko-')
axis([min(d_hkl(peaks2use)) max(d_hkl(peaks2use)) -1e-3 1e-3])
xlabel('d_{hkl} [Angstrom]')
ylabel('pseudo-strain [-]')
grid on
legend('using initial prms', 'using optimized prms')

figure,
plot(d_hkl(peaks2use), pseudo_strain, 'ko--')
xlabel('d_{hkl} [Angstrom]')
ylabel('pseudo-strain (-)')
grid on
axis([min(d_hkl(peaks2use)) max(d_hkl(peaks2use)) -1e-3 1e-3])
% axis square tight