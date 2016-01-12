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
ms              = load('fcc.ms');
DetectedPeaks   = [1 2 3 4];
qsym            = CubSymmetries;
PointsPerFiber  = 100;

ws  = load('wscub5x');
[dh, eVals] = DiscreteHarmonics(ws.wscub.frmesh, 23);

odf.frmesh  = ws.wscub.frmesh;
% odf.field   = ones(odf.frmesh.numind,1);
% odf.field   = 5*rand(odf.frmesh.numind,1);
odf.field   = dh(:,4) + abs(min(dh(:,3))); 
% odf.field   = odf.field + abs(min(odf.field));
% odf.field   = odf.field./MeanValue(odf.field, odf.frmesh.l2ip);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONVERT TO mm3
AverageGrainVolume  = AverageGrainVolume / 1000 / 1000 / 1000;

GaugeLengthZUS  = IncSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZDS  = OutSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZ    = GaugeLengthZUS + GaugeLengthZDS;
GaugeVolume     = GaugeLengthZ * IncSlitSizeAzi;            %%% THIS ASSUMES A TETRAGONAL BOX
NumberOfGrains  = GaugeVolume / AverageGrainVolume;

disp(sprintf('IncSlit = %3.3f mm & OutSlit = %3.3f mm', IncSlitSizeRad, OutSlitSizeRad))
disp(sprintf('Gauge Volume                          = %3.3f mm^3', GaugeVolume))
disp(sprintf('Average grain volume                  = %3.0f um^3', AverageGrainVolume * 1000^3))
disp(sprintf('Number of grains in the gauge volume  = %3.0e grains', NumberOfGrains))

%%% SCATTERING VECTORS IN LAB SYSTEM
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

%%% CRYSTAL FIBERS
for i = 1:1:length(DetectedPeaks)
    fibH{i} = FiberOfPoint(qH, hkls(:,DetectedPeaks(i)), PointsPerFiber, qsym);
    fibV{i} = FiberOfPoint(qV, hkls(:,DetectedPeaks(i)), PointsPerFiber, qsym);
    
    figure,
    PlotFRPerimeter('cubic');
    hold on
    plot3(fibH{i}(1,:), fibH{i}(2,:), fibH{i}(3,:), 'r.')
    plot3(fibV{i}(1,:), fibV{i}(2,:), fibV{i}(3,:), 'b.')
    axis square off 
    view([-45 20])
end

[MeanValue(odf.field, odf.frmesh.l2ip) SumValue(odf.field, odf.frmesh.l2ip)];
PlotFR(odf.frmesh, odf.field)

% pdf = odf.field./MeanValue(odf.field, odf.frmesh.l2ip);     %%% PROBABILITY DENSITY FUNCTION FROM ODF
% [MeanValue(pdf, odf.frmesh.l2ip) SumValue(pdf, odf.frmesh.l2ip)]
% PlotFR(odf.frmesh, pdf)

pdf = odf.field./SumValue(odf.field, odf.frmesh.l2ip);     %%% PROBABILITY DENSITY FUNCTION FROM ODF
[MeanValue(pdf, odf.frmesh.l2ip) SumValue(pdf, odf.frmesh.l2ip)];
PlotFR(odf.frmesh, pdf)

for i = 1:1:length(DetectedPeaks)    
    odfpf   = BuildOdfPfMatrix(hkls(:,DetectedPeaks(i)), ...
        odf.frmesh, odf.frmesh.symmetries, ...
        [qH qV], PointsPerFiber, 0, 500);
    
    nGrains = odfpf*pdf*NumberOfGrains/ms(DetectedPeaks(i));
    disp(sprintf('Number of grains interrogated by q(%2.3f, %2.3f, %2.3f) || c{%d %d %d} = %3.2e grains', ...
        qH(1), qH(2), qH(3),  hkls(1,DetectedPeaks(i)), hkls(2,DetectedPeaks(i)), hkls(3,DetectedPeaks(i)), nGrains(1)))
    disp(sprintf('Number of grains interrogated by q(%2.3f, %2.3f, %2.3f) || c{%d %d %d} = %3.2e grains', ...
        qV(1), qV(2), qV(3),  hkls(1,DetectedPeaks(i)), hkls(2,DetectedPeaks(i)), hkls(3,DetectedPeaks(i)), nGrains(2)))
end