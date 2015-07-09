clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOA  = 6.97005;
ChToEnergyConversion  = [0.0925699 -0.0754175];
MeasurementPlane    = 'h';

%%% NUMBER OF POINTS IN X,Y,Z (GET FROM PYTHON FILE)
NX  = 1; NY  = 1; NZ  = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCC Fe
latticeParms    = 2.868 ;                        % IN Angstrom
hkls            = load('bcc.hkls');
d_hkl           = PlaneSpacings(latticeParms, 'cubic', hkls');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAR FILE DESIGNATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname_pypar     = './strain-examples/';
fname_pypar     = 'mach_feb15_TOA_7_Lap7_1.pypar';
pfname_pypar    = fullfile(pname_pypar, fname_pypar);
pardata         = ReadPythonParFile(pfname_pypar, 'Version', 'feb15');
% pardata         = ReadPythonParFile(pfname_pypar, 'Version', 'mpe_standard');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH WHERE DATA FILES LIVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname_data  = './strain-examples/Lap7/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_hkl  = 2.*d_hkl*sind(TOA/2);
E_hkl       = Angstrom2keV(lambda_hkl);

E_grid  = 1:1:2048;
E_grid  = ChToEnergyConversion(1)*E_grid + ChToEnergyConversion(2);

numDV   = length(pardata.day);
data    = zeros(2048,1);
data_stack  = zeros(2048,numDV);
for i = 1:1:numDV
    if strcmpi(MeasurementPlane, 'h')
        fname   = pardata.HorzFileName{i}(1:end-1);
    elseif strcmpi(MeasurementPlane, 'v')
        fname   = pardata.VertFileName{i}(1:end-1);
    end
    
    pfname  = fullfile(pname_data, fname);
    datai	= load(pfname);
    data    = data + datai;
    data_stack(:,i)   = datai;
end

if strcmpi(MeasurementPlane, 'h')
    titlestr    = 'data from horizontal detector';
elseif strcmpi(MeasurementPlane, 'v')
    titlestr    = 'data from vertical detector';
end

figure(1)
subplot(2,1,1)
imagesc(data_stack)
YLabelTick  = str2num(get(gca, 'YTickLabel'));
YLabelTick  = round(ChToEnergyConversion(1)*YLabelTick + ChToEnergyConversion(2));
YLabelTick  = num2str(YLabelTick);
set(gca, 'YTickLabel', YLabelTick)
grid on
xlabel('data number')
ylabel('E (keV)')
title(titlestr)

subplot(2,1,2)
plot(E_grid, data, 'k.')
hold on
plot(E_hkl, ones(length(E_hkl), 1), 'r^')
grid on
axis([min(E_grid) max(E_grid) min(data) max(data)])
xlabel('Energy (keV)')
ylabel('counts')
title(titlestr)
