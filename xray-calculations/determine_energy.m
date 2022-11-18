clear all
close all
clc

% addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

E_grid_keV  = 40:0.1:130;
E_grid_MeV  = E_grid_keV./1000;

[MAC_Ni, ~, density_Ni] = PhotonAttenuation('Ni', E_grid_MeV);
[MAC_Sn, ~, density_Sn] = PhotonAttenuation('Sn', E_grid_MeV);

sample_thickness_mm = 1.2;
sample_thickness_cm = sample_thickness_mm./10;

transmission_pct_Sn = exp(-MAC_Sn*density_Sn*sample_thickness_cm);
transmission_pct_Ni = exp(-MAC_Ni*density_Ni*sample_thickness_cm);

figure,
semilogy(E_grid_keV, transmission_pct_Ni)
hold on
semilogy(E_grid_keV, transmission_pct_Sn)
grid on
axis([60 90 0.1 1])
xlabel('E [keV]')
ylabel('transmission [%]')
legend('Ni', 'Sn')

% hold on
% semilogy(edd_bm_fluxdata_raw(:,1), edd_bm_fluxdata_raw(:,2).*edd_tr(:,1))
% semilogy(edd_bm_fluxdata_raw(:,1), edd_bm_fluxdata_raw(:,2).*edd_tr(:,2))
% semilogy(edd_bm_fluxdata_raw(:,1), edd_bm_fluxdata_raw(:,2).*edd_tr(:,3))
% legend('raw', '1x ss-filter', '2x ss-filter', '3x ss-filter')
% grid on
% axis([20 250 1 1e15])
% xlabel('E [keV]')
% ylabel('flux [ph/s/0.1%bw]')


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EDD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edd_bm_fluxdata_raw = load('bm_flux.data');

E_grid_edd_bm_flux_keV  = edd_bm_fluxdata_raw(:,1);
E_grid_edd_bm_flux_MeV  = E_grid_edd_bm_flux_keV./1000;
[MAC_Sn, ~, density_Sn]       = PhotonAttenuation('Fe', E_grid_edd_bm_flux_MeV);

filter_thickness_mm(1)  = 0.7;
filter_thickness_mm(2)  = 1.4;
filter_thickness_mm(3)  = 2.1;
filter_thickness_cm     = filter_thickness_mm./10;

edd_tr(:,1) = exp(-MAC_Sn*density_Sn*filter_thickness_cm(1));
edd_tr(:,2) = exp(-MAC_Sn*density_Sn*filter_thickness_cm(2));
edd_tr(:,3) = exp(-MAC_Sn*density_Sn*filter_thickness_cm(3));

figure,
semilogy(edd_bm_fluxdata_raw(:,1), edd_bm_fluxdata_raw(:,2))
hold on
semilogy(edd_bm_fluxdata_raw(:,1), edd_bm_fluxdata_raw(:,2).*edd_tr(:,1))
semilogy(edd_bm_fluxdata_raw(:,1), edd_bm_fluxdata_raw(:,2).*edd_tr(:,2))
semilogy(edd_bm_fluxdata_raw(:,1), edd_bm_fluxdata_raw(:,2).*edd_tr(:,3))
legend('raw', '1x ss-filter', '2x ss-filter', '3x ss-filter')
grid on
axis([20 250 1 1e15])
xlabel('E [keV]')
ylabel('flux [ph/s/0.1%bw]')

for iii = 1:1:3
    edd_tr_ratio(:,iii) = edd_tr(:,1)./edd_tr(:,iii);
end

% edd_flux(1,:)       = interp1(fluxdata(:,1), fluxdata(:,2), Efit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DETECTOR CALIBRATION INFORMATION
%%% det_id  = 1;
ChToEnergyConversion(1,:)   =  [0.025622209096602   0.186706774011202]; % det1
TOA(1)                      = 6.470742865420140;

%%% det_id  = 2;
ChToEnergyConversion(2,:)   = [0.026521679794167   0.011642807382541]; % det2
TOA(2)                      = 6.875344433681863;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xgrid       = 1:1:8192;
Egrid_edd(:,1)  = polyval(ChToEnergyConversion(1,:), xgrid);
Egrid_edd(:,2)  = polyval(ChToEnergyConversion(2,:), xgrid);

lambda_grid_edd(:,1)    = keV2Angstrom(Egrid_edd(:,1));
lambda_grid_edd(:,2)    = keV2Angstrom(Egrid_edd(:,2));

dspacing_grid_edd(:,1)  = lambda_grid_edd(:,1)./(2*sind(TOA(1)./2));
dspacing_grid_edd(:,2)  = lambda_grid_edd(:,2)./(2*sind(TOA(2)./2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% METADATA
% pfname_metadata_edd = '/home/s1b/__eval/projects_parkjs/edd_6bma_AFRL_Elmer_oct21/exposure_record_full.par';
pfname_metadata_edd = 'D:\w\edd_6bma_AFRL_Elmer_oct21\exposure_record_full.par';
metadata_edd        = ReadSpecParFile(pfname_metadata_edd, 'Version', 'mpe_standard');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FeCu5050
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proot_edd   = '/home/s1b/__eval/projects_parkjs/edd_6bma_AFRL_Elmer_oct21';
proot_edd   = 'D:\w\edd_6bma_AFRL_Elmer_oct21';

froot_edd   = 'Cu5050_100x100um_1filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_5050(1) = readedd5_6bm(pfroot);

froot_edd   = 'Cu5050_100x100um_2filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_5050(2) = readedd5_6bm(pfroot);

froot_edd   = 'Cu5050_100x100um_3filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_5050(3) = readedd5_6bm(pfroot);

for iii = 1:1:length(edd_5050)
    figure(iii)
    for detid = 1:1:2
        semilogy(dspacing_grid_edd(:,detid), edd_5050(iii).data{detid}', '-')
        hold on
    end
    xlabel('d-spacing [Angstrom]')
    ylabel('intensity [cts]')
    axis([0.5 2.3 0 1e3])
end

%%% intensity correction for ss filter
for iii = 1:1:length(edd_5050)
    for detid = 1:1:2
        
        intensity_multiplier_ss_filter  = interp1(E_grid_edd_bm_flux_keV, edd_tr_ratio(:,iii), Egrid_edd(:,detid));
        
        edd_5050(iii).intensity_multiplier_ss_filter{detid} = intensity_multiplier_ss_filter;
        edd_5050(iii).ss_filter_corrected_data{detid}       = edd_5050(iii).data{detid}'.*intensity_multiplier_ss_filter;
        
        figure(detid+10)
        % semilogy(dspacing_grid(:,detid), edd_5050(iii).data{detid}'.*intensity_multiplier, '-')
        plot(dspacing_grid_edd(:,detid), edd_5050(iii).ss_filter_corrected_data{detid}, '-')
        hold on
        title(sprintf('intensity corrected for ss filter det-%d', detid))
        xlabel('d-spacing [Angstrom]')
        ylabel('intensity [cts]')
        axis([0.5 2.3 0 3e3])
    end
end

%%% FeCu6733
froot_edd   = 'Cu6733_100x100um_1filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_6733(1) = readedd5_6bm(pfroot);

froot_edd   = 'Cu6733_100x100um_2filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_6733(2) = readedd5_6bm(pfroot);

froot_edd   = 'Cu6733_100x100um_3filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_6733(3) = readedd5_6bm(pfroot);

for iii = 1:1:length(edd_6733)
    for detid = 1:1:2
        
        intensity_multiplier_ss_filter  = interp1(E_grid_edd_bm_flux_keV, edd_tr_ratio(:,iii), Egrid_edd(:,detid));
        
        edd_6733(iii).intensity_multiplier_ss_filter{detid} = intensity_multiplier_ss_filter;
        edd_6733(iii).ss_filter_corrected_data{detid}       = edd_6733(iii).data{detid}'.*intensity_multiplier_ss_filter;
        
        figure(detid+20)
        % semilogy(dspacing_grid(:,detid), edd_5050(iii).data{detid}'.*intensity_multiplier, '-')
        plot(dspacing_grid_edd(:,detid), edd_6733(iii).ss_filter_corrected_data{detid}, '-')
        hold on
        title(sprintf('intensity corrected for ss filter det-%d', detid))
        xlabel('d-spacing [Angstrom]')
        ylabel('intensity [cts]')
        axis([0.5 2.3 0 3e3])
    end
end

%%% CuRD
froot_edd   = 'CuRD_100x100um_1filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
tempedd     = readedd5_6bm(pfroot);
edd_CuRD(1) = tempedd(2);

froot_edd   = 'CuRD_100x100um_2filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_CuRD(2) = readedd5_6bm(pfroot);

froot_edd   = 'CuRD_100x100um_3filter';
pfroot      = fullfile(proot_edd, froot_edd, sprintf('%s.xy', froot_edd));
edd_CuRD(3) = readedd5_6bm(pfroot);

for iii = 1:1:length(edd_CuRD)
    for detid = 1:1:2
        
        intensity_multiplier_ss_filter  = interp1(E_grid_edd_bm_flux_keV, edd_tr_ratio(:,iii), Egrid_edd(:,detid));
        
        edd_CuRD(iii).intensity_multiplier_ss_filter{detid} = intensity_multiplier_ss_filter;
        edd_CuRD(iii).ss_filter_corrected_data{detid}       = edd_CuRD(iii).data{detid}'.*intensity_multiplier_ss_filter;
        
        figure(detid+30)
        % semilogy(dspacing_grid(:,detid), edd_5050(iii).data{detid}'.*intensity_multiplier, '-')
        plot(dspacing_grid_edd(:,detid), edd_CuRD(iii).ss_filter_corrected_data{detid}, '-')
        hold on
        title(sprintf('intensity corrected for ss filter det-%d', detid))
        xlabel('d-spacing [Angstrom]')
        ylabel('intensity [cts]')
        axis([0.5 2.3 0 1.6e4])
    end
end
