clear all
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MaterialName    = 'Al';                         % FCC Al 
% latticeParms    = 4.050;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% [a, b, c]       = StructureFactor(MaterialName);
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Ce';                         % BCC Fe
% latticeParms    = 5.4114 ;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Fe';                         % BCC Fe
% latticeParms    = 2.87 ;                        % IN Angstrom
% hkls            = load('bcc.hkls');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Fe';                         % FCC Fe
% latticeParms    = 3.515;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Ni';                         % FCC Ni
% latticeParms    = 3.520;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'C';                         % diamond
% latticeParms    = 3.5668;                        % IN Angstrom
% hkls            = load('diamondcubic.hkls');
% [a, b, c]       = StructureFactor(MaterialName);
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

MaterialName    = 'Ti';                         % diamond
latticeParms    = [2.95111 4.68433];                        % IN Angstrom
hkls            = load('hcp.hkls');
[a, b, c]       = StructureFactor(MaterialName);
d_hkls          = PlaneSpacings(latticeParms, 'hexagonal', hkls');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SampleThickness = 5.0;                        % IN cm
TakeOffAngle    = 2:0.1:7.0;                  % IN deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BeamLineFlux    = load('bm_flux.data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numhkls         = size(hkls,1);

%%% STRUCTURE FACTOR CALCULATION
f   = c;
x   = 0.5./d_hkls;
for i = 1:1:4
    f   = f + a(i)*exp(-b(i)*x.*x);
end

%%% FOR FCC
F   = 4*f';
%%% FOR BCC
% F   = 2*f';
%%% FOR SIMPLE CUBIC
% F   = f';

F2	= F.^2;
F2_normalized   = F2./max(F2);
F2_normalized   = ones(numhkls, 1);

%%% INTENSITY CALCULATION
for j = 1:1:length(TakeOffAngle)
    wavelength      = 2*d_hkls.*sind(TakeOffAngle(j)/2);
    Energy          = Angstrom2keV(wavelength);         % IN keV
    
    [mu, ~, rho]    = PhotonAttenuation(MaterialName, Energy'./1000);       % cm^2/g
    
    PercentTransmissionDS(:,j)    = exp(-mu*rho*SampleThickness);
    PercentTransmissionUS(:,j)    = exp(-mu*rho*(SampleThickness/cosd(TakeOffAngle(j))));
    
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
    PhotonTransAtNormalIncDS(:,j)       = PercentTransmissionDS(:,j).*F2_normalized.*Flux;
    PhotonTransAtNormalIncUS(:,j)       = PercentTransmissionUS(:,j).*F2_normalized.*Flux;
end

figure,
set(gcf, 'Position', [1000 50 900 950])
subplot(2,3,1)
semilogy(BeamLineFlux(:,1), BeamLineFlux(:,2), 'k.')
xlabel('Energy (keV)')
ylabel('Number of photons')
title('BM Flux (photons / s / 0.1% BW)')
grid on

subplot(2,3,2)
plot(1:1:numhkls, TransEnergy(:,1), 'b.')
hold on
plot(1:1:numhkls, TransEnergy(:,end), 'r.')
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Energy (keV)')
title('Diffraction energy - First and the last TOA only')
grid on

subplot(2,3,3)
plot(1:1:numhkls, F, 'k.')
xlabel('hkl id')
ylabel('Structure Factor (-)')
grid on

subplot(2,3,4)
plot(1:1:numhkls, PercentTransmissionDS(:,1)*100, 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5)
hold on
plot(1:1:numhkls, PercentTransmissionDS(:,end)*100, 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
plot(1:1:numhkls, PercentTransmissionUS(:,1)*100, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
plot(1:1:numhkls, PercentTransmissionUS(:,end)*100, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
% axis([1 numhkls 0 0.001])
legend(sprintf('DS 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('DS 2th = %2.1f deg', TakeOffAngle(end)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Percent transmission')
grid on

subplot(2,3,5)
semilogy(1:1:numhkls, PhotonTransAtNormalIncDS(:,1), 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5)
hold on
semilogy(1:1:numhkls, PhotonTransAtNormalIncDS(:,end), 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
semilogy(1:1:numhkls, PhotonTransAtNormalIncUS(:,1), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
semilogy(1:1:numhkls, PhotonTransAtNormalIncUS(:,end), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
% axis([1 numhkls 1e-7 1e14])
legend(sprintf('DS 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('DS 2th = %2.1f deg', TakeOffAngle(end)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Number of photons transmitted')
grid on
% axis([1 numhkls 10^8 10^14])
return
figure,
set(gcf, 'Position', [1000 50 900 950])
subplot(2,3,1)
semilogy(BeamLineFlux(:,1), BeamLineFlux(:,2), 'k.')
xlabel('Energy (keV)')
ylabel('Number of photons')
title('BM Flux (photons / s / 0.1% BW)')
grid on

subplot(2,3,2)
plot(d_hkls, TransEnergy(:,1), 'b.')
hold on
plot(d_hkls, TransEnergy(:,end), 'r.')
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('d-spacing (Angstrom)')
ylabel('Energy (keV)')
title('Diffraction energy - First and the last TOA only')
grid on

subplot(2,3,3)
plot(d_hkls, F(:,1), 'k.')
xlabel('d-spacing (Angstrom)')
ylabel('Structure Factor (-)')
grid on

subplot(2,3,4)
plot(d_hkls, PercentTransmissionDS(:,1)*100, 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5)
hold on
plot(d_hkls, PercentTransmissionDS(:,end)*100, 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
plot(d_hkls, PercentTransmissionUS(:,1)*100, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
plot(d_hkls, PercentTransmissionUS(:,end)*100, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
% axis([1 numhkls 0 0.001])
legend(sprintf('DS 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('DS 2th = %2.1f deg', TakeOffAngle(end)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(end)), 'Location', 'Best')
xlabel('d-spacing (Angstrom)')
ylabel('Percent transmission')
grid on

subplot(2,3,5)
semilogy(d_hkls, PhotonTransAtNormalIncDS(:,1), 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5)
hold on
semilogy(d_hkls, PhotonTransAtNormalIncDS(:,end), 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
semilogy(d_hkls, PhotonTransAtNormalIncUS(:,1), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5)
semilogy(d_hkls, PhotonTransAtNormalIncUS(:,end), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
% axis([1 numhkls 1e-7 1e14])
legend(sprintf('DS 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('DS 2th = %2.1f deg', TakeOffAngle(end)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(1)), ...
    sprintf('US 2th = %2.1f deg', TakeOffAngle(end)), 'Location', 'Best')
xlabel('d-spacing (Angstrom)')
ylabel('Number of photons transmitted')
grid on