clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TakeOffAngle    = 3:0.5:7.0;                    % IN deg
IncSlitSizeRad  = 0.2;                          % IN mm
OutSlitSizeRad  = 0.2;                          % IN mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GaugeLengthZUS  = IncSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZDS  = OutSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZ    = GaugeLengthZUS + GaugeLengthZDS;

figure(2)
plot(TakeOffAngle, GaugeLengthZ, 'ko-')
axis tight
title(['IncSlit = ', num2str(IncSlitSizeRad), ' mm & OutSlit = ', num2str(OutSlitSizeRad), ' mm'])
xlabel('take off angle (deg)')
ylabel('gauge length in z (mm)')
grid on
