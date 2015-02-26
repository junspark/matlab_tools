clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOA  = 6.97005;
ChToEnergyConversion  = [0.0925699 -0.0754175];
MeasurementPlane    = 'h';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCC Fe
latticeParms    = 2.868 ;                        % IN Angstrom
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
% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XS  = pardata.samX;
YS  = pardata.samY;
ZS  = pardata.samZ;

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
end
return
for i = 1:1:length(pkid_fit)
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
end
