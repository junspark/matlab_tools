clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MaterialName    = 'Al';                         % FCC Al 
% latticeParms    = 4.050;                        % IN Angstrom
% hkls            = load('fcc.hkls');

MaterialName    = 'Fe';                         % BCC Fe
latticeParms    = 2.87 ;                        % IN Angstrom
hkls            = load('bcc.hkls');

% MaterialName    = 'Fe';                         % FCC Fe
% latticeParms    = 3.515;                        % IN Angstrom
% hkls            = load('fcc.hkls');

% MaterialName    = 'Ni';                         % FCC Ni
% latticeParms    = 3.520;                        % IN Angstrom
% hkls            = load('fcc.hkls');

% MaterialName    = 'C';                         % diamond
% latticeParms    = 3.5668;                        % IN Angstrom
% hkls            = load('diamondcubic.hkls');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SampleThickness = 30;                          % IN cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IncSlitSizeRad  = 0.2;                          % IN mm
OutSlitSizeRad  = 0.2;                          % IN mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numhkls         = size(hkls,1);
d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

BeamLineFlux    = load('bm_flux.data');
TakeOffAngle    = 5:1:10;                    % IN deg

GaugeLengthZUS  = IncSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZDS  = OutSlitSizeRad*cosd(TakeOffAngle./2)./sind(TakeOffAngle);
GaugeLengthZ    = GaugeLengthZUS + GaugeLengthZDS

for j = 1:1:length(TakeOffAngle)
    wavelength      = 2*d_hkls.*sind(TakeOffAngle(j)/2);
    Energy          = Angstrom2keV(wavelength);         % IN keV
    
    [mu, ~, rho]    = PhotonAttenuation(MaterialName, Energy'./1000);       % cm^2/g
    
    PercentTransmission(:,j)    = exp(-mu*rho*SampleThickness);
    
    Flux    = zeros(length(Energy),1);
    for i = 1:1:length(Energy)
        idx = find(Energy(i) == BeamLineFlux(:,1));
        if isempty(idx)
            idx1    = find(Energy(i) > BeamLineFlux(:,1));
            idx2    = find(Energy(i) < BeamLineFlux(:,1));
            if isempty(idx1) || isempty(idx2)
                Flux(i,1)   = nan;
            else
                idx1    = idx1(end);
                idx2    = idx2(1);
                
                Flux(i,1)   = BeamLineFlux(idx1,2) - (BeamLineFlux(idx1,2) - BeamLineFlux(idx2,2))/(BeamLineFlux(idx1,1) - BeamLineFlux(idx2,1)) * (BeamLineFlux(idx1,1) - Energy(i));
            end
        else
            Flux(i,1)   = BeamLineFlux(idx,2);
        end
    end
    
    TransEnergy(:,j)                    = Energy';
    PhotonTransAtNormalIncidence(:,j)   = PercentTransmission(:,j).*Flux;
end

figure,
set(gcf, 'Position', [1007 33 902 977])
subplot(2,2,1)
semilogy(BeamLineFlux(:,1), BeamLineFlux(:,2), 'k.')
xlabel('Energy (keV)')
ylabel('Number of photons')
title('BM Flux (photons / s / 0.1% BW)')
grid on

subplot(2,2,2)
plot(1:1:numhkls, TransEnergy(:,1), 'b.')
hold on
plot(1:1:numhkls, TransEnergy(:,end), 'r.')
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Energy (keV)')
title('Diffraction energy - First and the last TOA only')
grid on

subplot(2,2,3)
plot(1:1:numhkls, PercentTransmission(:,1)*100, 'b.')
hold on
plot(1:1:numhkls, PercentTransmission(:,end)*100, 'r.')
% axis([1 numhkls 0 0.001])
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Percent transmission')
grid on

subplot(2,2,4)
semilogy(1:1:numhkls, PhotonTransAtNormalIncidence(:,1), 'b.')
hold on
semilogy(1:1:numhkls, PhotonTransAtNormalIncidence(:,end), 'r.')
% axis([1 numhkls 1e-7 1e14])
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Number of photons transmitted')
grid on

figure,
plot(TakeOffAngle, GaugeLengthZ, 'ko')
axis tight
title(['IncSlit = ', num2str(IncSlitSizeRad), ' mm & OutSlit = ', num2str(IncSlitSizeRad), ' mm'])
xlabel('take off angle (deg)')
ylabel('gauge length in z (mm)')
grid on