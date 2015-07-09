clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTAL GEOMETRY
% ASSUME FOR NOW A POINT DETECTOR (NO AZIMUTHAL ANGLE)
% THIS IS THE LOWER BOUND FOR THE NUMBER OF GRAINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TakeOffAngle    = 5;                            % IN deg
IncSlitSizeRad  = 0.2;                          % IN mm
OutSlitSizeRad  = 0.2;                          % IN mm

IncSlitSizeAzi  = 0.2;                          % IN mm
OutSlitSizeAzi  = 0.2;                          % IN mm

OutSlitDsam     = 100;                          % IN mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL & SAMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AverageGrainVolume  = 100;                      % IN um^3
hkls            = load('fcc.hkls')';
DetectedPeaks   = [1 2 3 4];
qsym            = CubSymmetries;

PointsPerFiber  = 1000;

ws  = load('wscub4x');
odf.frmesh      = ws.wscub4x.frmesh;
odf.field       = ones(odf.frmesh.numind,1);

% CONVERT TO mm3
AverageGrainVolume  = AverageGrainVolume / 1000 / 1000 / 1000;

GaugeLengthZUS  = IncSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZDS  = OutSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZ    = GaugeLengthZUS + GaugeLengthZDS;
GaugeVolume     = GaugeLengthZ * IncSlitSizeAzi;
NumberOfGrains  = GaugeVolume / AverageGrainVolume;

disp(sprintf('IncSlit = %3.3f mm & OutSlit = %3.3f mm', IncSlitSizeRad, OutSlitSizeRad))
disp(sprintf('Gauge Volume                          = %3.3f mm^3', GaugeVolume))
disp(sprintf('Average grain volume                  = %3.0f um^3', AverageGrainVolume * 1000^3))
disp(sprintf('Number of grains in the gauge volume  = %3.0e grains', NumberOfGrains))

%%% SCATTERING VECTOR IN LAB SYSTEM
Ry  = [ ...
    cosd(TakeOffAngle/2) 0 sind(TakeOffAngle/2); ...
    0 1 0; ...
    -sind(TakeOffAngle/2) 0 cosd(TakeOffAngle/2); ...
    ];

Rx  = [ ...
    1 0 0; ...
    0 cosd(-TakeOffAngle/2) -sind(-TakeOffAngle/2); ...
    0 sind(-TakeOffAngle/2) cosd(-TakeOffAngle/2); ...
    ];

qH  = Ry*[1 0 0]';
qV  = Rx*[0 1 0]';

for i = 1:1:length(DetectedPeaks)
    fib_H   = FiberOfPoint(qH, hkls(:,DetectedPeaks(i)), PointsPerFiber, qsym);
    [fele, fcrd]    = FiberCoordinates(fib_H, odf.frmesh);
    values          = EvalMeshFunc(odf.frmesh, odf.field, fele, fcrd);
    sum(values)
    retutn
    figure(i)
    PlotFRPerimeter('cubic')
    axis off equal
    plot3(fib_H(1,:), fib_H(2,:), fib_H(3,:), 'k.')
end
return
qV  = [];

figure(2)
plot(TakeOffAngle, GaugeLengthZ, 'ko-')
axis tight
title(['IncSlit = ', num2str(IncSlitSizeRad), ' mm & OutSlit = ', num2str(IncSlitSizeRad), ' mm'])
xlabel('take off angle (deg)')
ylabel('gauge length in z (mm)')
grid on