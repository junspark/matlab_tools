%%% CAN BE USED TO LOOK AT MOTOR_CAM STACK FOR ALIGNMENT

clear all
close all
clc

%%% BASED ON SIMULATION WE SHOULD HAVE 4.17 mm ALONG Z
% vertical center   : 1.619024
% vertical fwhm     : 3.188952 
% horizontal center : 1.612883
% horizontal fwhm   : 4.181930
% 50th point to use for TOA calculation

proot   = '/home/s1b/__eval/projects_parkjs/edd_6bm_2017-1/startup_mar17';

format long

%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%
pname     = 'CeO2_hv_gv_20170314';

pt_ini  = -3.3;
pt_fin  = 6.7;
npts    = 101;
%%%%%

pts_grid    = linspace(pt_ini, pt_fin, npts);
froot       = pname;

data_v  = zeros(8192,npts);
data_h  = zeros(8192,npts);

data_v_sum = zeros(npts,1);
data_h_sum = zeros(npts,1);

roi_v   = 1:8192;
roi_h   = roi_v;

for i = 1:1:npts
    fname   = sprintf('%s-%03d-hv.xy', froot, i);
    pfname  = fullfile(proot, pname, fname);
    while (exist(pfname) ~= 2)
        disp(sprintf('waiting for %s', fname))
        pause(2)
    end
    disp(fname);
    
    data    = load(pfname);
    data_v(:,i) = data(:,1);
    data_h(:,i) = data(:,2);
    
    data_v_sum(i)   = sum(data(roi_v,1));
    data_h_sum(i)   = sum(data(roi_h,2));
end

pkfitting_pars.pfunc_type   = 'pseudovoigt';
pkfitting_pars.pbkg_order   = 2;

%%% V FIT
p0  = [ max(data_v_sum); std(pts_grid);   0.5; mean(pts_grid); 0; mean(data_v_sum); ];
pLB = [          0;   0;    0;    min(pts_grid); -inf; -inf; ];
pUB = [        inf; inf;    1;    max(pts_grid);  inf;  inf; ];
pkfitting_pars.xdata    = pts_grid;

[pv, rnv, ~, efv]    = lsqcurvefit(@pfunc_switch, ...
    p0, pkfitting_pars, data_v_sum, pLB, pUB);
vfit0   = pfunc_switch(p0, pkfitting_pars);
vfit    = pfunc_switch(pv, pkfitting_pars);
disp(sprintf('vertical center   : %f', pv(4)));
disp(sprintf('vertical fwhm     : %f', pv(2)));
 
%%% H FIT
p0  = [ max(data_h_sum); std(pts_grid);   0.5;   mean(pts_grid); 0; mean(data_h_sum); ];
pLB = [          0;   0;    0;    min(pts_grid); -inf; -inf; ];
pUB = [        inf; inf;    1;    max(pts_grid);  inf;  inf; ];
pkfitting_pars.xdata    = pts_grid;

[ph, rnh, ~, efh]    = lsqcurvefit(@pfunc_switch, ...
    p0, pkfitting_pars, data_h_sum, pLB, pUB);
hfit0   = pfunc_switch(p0, pkfitting_pars);
hfit    = pfunc_switch(ph, pkfitting_pars);
disp(sprintf('horizontal center : %f', ph(4)));
disp(sprintf('horizontal fwhm   : %f', ph(2)));

figure,
subplot(2,1,1);
imagesc(log(data_v'))
title(sprintf('%s-v', pname), 'interpreter', 'none')
xlabel('chan number')
ylabel('scan number')
colorbar vert
colormap jet
% caxis([0 100])

subplot(2,1,2);
imagesc((data_h'))
title(sprintf('%s-h', pname), 'interpreter', 'none')
xlabel('chan number')
ylabel('scan number')
colorbar vert
colormap jet
% caxis([0 100])

figure,
subplot(1,2,1);
plot(pts_grid, data_v_sum, 'b+')
hold on
plot(pts_grid, vfit0, 'r:')
plot(pts_grid, vfit, 'g-')
title(sprintf('%s-v', pname), 'interpreter', 'none')
xlabel('grid (phys units)')
ylabel('cts')

subplot(1,2,2);
plot(pts_grid, data_h_sum, 'b+')
hold on
plot(pts_grid, hfit0, 'r:')
plot(pts_grid, hfit, 'g-')
title(sprintf('%s-h', pname), 'interpreter', 'none')
xlabel('grid (phys units)')
ylabel('cts')
