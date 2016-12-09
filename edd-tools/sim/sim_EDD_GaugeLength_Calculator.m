clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TakeOffAngle    = 2:0.5:5.0;                    % IN deg
IncSlitSizeRad  = 0.1;                          % IN mm
OutSlitSizeRad  = 0.1;                          % IN mm
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
