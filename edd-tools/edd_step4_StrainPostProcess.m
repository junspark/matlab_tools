clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOA     =  6.970451365881403;
ChToEnergyConversion    = [0.092570164886107  -0.077311622210394];
MeasurementPlane        = 'h';

%%% NUMBER OF POINTS IN X,Y,Z (GET FROM PYTHON FILE)
NX  = 15; NY  = 1; NZ  = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCC Fe
latticeParms    = 2.8698 ;                        % IN Angstrom
hkls            = load('bcc.hkls');
d_hkl           = PlaneSpacings(latticeParms, 'cubic', hkls');
pkid_fit        = 4:8;
sqrt_hkls       = sqrt(sum(hkls(pkid_fit, :).*hkls(pkid_fit, :),2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAR FILE DESIGNATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname_pypar     = './strain-examples/';
fname_pypar     = 'mach_feb15_TOA_7_Lap7.pypar';
pfname_pypar    = fullfile(pname_pypar, fname_pypar);
pardata         = ReadPythonParFile(pfname_pypar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH WHERE DATA FILES LIVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname_data  = './strain-examples/Lap7/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE NAME FOR STRAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname_strain    = 'mach_feb15_strain_h';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XS  = pardata.samX;
YS  = pardata.samY;
ZS  = pardata.samZ;

xs  = reshape(XS', NX, NY, NZ);
ys  = reshape(YS', NX, NY, NZ);
zs  = reshape(ZS', NX, NY, NZ);

lambda_hkl  = 2.*d_hkl*sind(TOA/2);
E_hkl       = Angstrom2keV(lambda_hkl);

E_grid  = 1:1:2048;
E_grid  = ChToEnergyConversion(1)*E_grid + ChToEnergyConversion(2);

numDV   = length(pardata.day);
for i = 1:1:numDV
    if strcmpi(MeasurementPlane, 'h')
        fname   = pardata.HorzFileName{i}(1:end-1);
    elseif strcmpi(MeasurementPlane, 'v')
        fname   = pardata.VertFileName{i}(1:end-1);
    end
    
    fname_fit   = [fname, '.fit.mat'];
    pfname_fit  = fullfile(pname_data, fname_fit);
    
    fit_data    = load(pfname_fit);
    
    Efit(i,:)   = fit_data.Efit;
    Afit(i,:)   = fit_data.Afit;
    gfit(i,:)   = fit_data.gfit;
    nfit(i,:)   = fit_data.nfit;
    rn(i,:)     = fit_data.rn;
    
    Rwp(i,:)    = fit_data.Rwp;
    Rp(i,:)     = fit_data.Rp;
    Re(i,:)     = fit_data.Re; 
end

for i = 1:1:length(pkid_fit)
    strain  = E_hkl(pkid_fit(i))./Efit(:,i) - 1;
    
    figure(100 + i)
    subplot(2,3,1)
    scatter3(XS, YS, ZS, 50, Afit(:,i), 'filled')
    colorbar vert
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    title('peak amplitude (cts)')
    
    subplot(2,3,2)
    scatter3(XS, YS, ZS, 50, gfit(:,i), 'filled')
    colorbar vert
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    title('peak width (keV)')
    
    subplot(2,3,3)
    scatter3(XS, YS, ZS, 50, nfit(:,i), 'filled')
    colorbar vert
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    title('mix parameter')
    
    subplot(2,3,4)
    scatter3(XS, YS, ZS, 50, Efit(:,i), 'filled')
    colorbar vert
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    title('E fit (keV)')
    
    subplot(2,3,5)
    scatter3(XS, YS, ZS, 50, Rwp(:,i), 'filled')
    colorbar vert
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    title('peak fit weigted residual')
    
    subplot(2,3,6)
    scatter3(XS, YS, ZS, 50, strain, 'filled')
    colorbar vert
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    title('strain')
    caxis([-3e-3 3e-3])
    
    Strain(:,i) = strain;
end

output  = [XS YS ZS Strain];
xlswrite(fname_strain, output)
