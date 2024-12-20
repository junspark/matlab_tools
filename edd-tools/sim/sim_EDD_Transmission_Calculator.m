clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MaterialName    = 'Al';                         % FCC Al 
% latticeParms    = 4.050;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% multiplicity    = load('fcc.ms');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Cu';                         % FCC Al 
% latticeParms    = 3.615;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% multiplicity    = load('fcc.ms');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Ce';                         % FCC Ce
% latticeParms    = 5.4114 ;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% multiplicity    = load('fcc.ms');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

MaterialName    = 'Fe';                         % BCC Fe
latticeParms    = 2.87 ;                        % IN Angstrom
hkls            = load('bcc.hkls');
multiplicity    = load('bcc.ms');
d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Fe';                         % FCC Fe
% latticeParms    = 3.515;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% multiplicity    = load('fcc.ms');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Ni';                         % FCC Ni
% latticeParms    = 3.520;                        % IN Angstrom
% hkls            = load('fcc.hkls');
% multiplicity    = load('fcc.ms');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'C';                         % diamond
% latticeParms    = 3.5668;                        % IN Angstrom
% hkls            = load('diamondcubic.hkls');
% multiplicity    = load('diamondcubic.ms');
% d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

% MaterialName    = 'Ti';                         % diamond
% latticeParms    = [2.95111 4.68433];                        % IN Angstrom
% hkls            = load('hcp.hkls');
% multiplicity    = load('hcp.ms');
% d_hkls          = PlaneSpacings(latticeParms, 'hexagonal', hkls');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SampleThickness = 2.5;                % IN cm
TakeOffAngle    = 3:6;                  % IN deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CONSTANTS
BeamLineFlux    = load('bm_flux.data');
GeThickness     = 1;    % ASSUME GE XSTAL THICKNESS 10 mm = 1 cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numhkls = size(hkls,1);
P       = 1; % POLARIZATION FACTOR (NO MONO; 1 FOR VERTICAL DETECTOR; ESSENTIALLY 1 AS BRAGG ANGLE SMALL)
L       = LorentzFactor(TakeOffAngle./2);
f       = AtomicScatteringFactor(MaterialName, d_hkls);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STRUCTURE FACTOR CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FOR FCC
F   = 4*f';
%%% FOR BCC
% F   = 2*f';
%%% FOR SIMPLE CUBIC
% F   = f';
%%% FOR HCP
% for i = 1:1:numhkls
%     h   = hkls(i,1);
%     k   = hkls(i,2);
%     l   = hkls(i,3);
%     h2k = h + 2*k;
%     
%     if mod(l,2) == 0
%         if mod(h2k,3) == 0
%             F(i,1)  = 2*f(i);
%         else
%             F(i,1)  = f(i);
%         end
%     else
%         if mod(h2k,3) == 0
%             F(i,1)  = 0;
%         else
%             F(i,1)  = sqrt(3)*f(i);
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTENSITY CALCULATION
F2	= F.^2;
F2_multiplicity     = F2.*multiplicity;
F2_multiplicity_LP  = (F2_multiplicity*P).*L;
for j = 1:1:length(TakeOffAngle)
    wavelength      = 2*d_hkls.*sind(TakeOffAngle(j)/2);
    Energy          = Angstrom2keV(wavelength);         % IN keV
    
    [mu, ~, rho]    = PhotonAttenuation(MaterialName, Energy'./1000);       % cm^2/g
    
    PercentTransmissionUS(:,j)    = exp(-mu*rho*(SampleThickness/cosd(TakeOffAngle(j))));
    PercentTransmissionDS(:,j)    = exp(-mu*rho*SampleThickness);
    
    [mu_Ge, ~, rho_Ge]          = PhotonAttenuation('Ge', Energy'./1000);       % cm^2/g
    PercentAbsorptionGe(:,j)    = 1 - exp(-mu_Ge*rho_Ge*GeThickness);
    
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
    
    F2_normalized(:,j)                  = F2_multiplicity_LP(:,j)./max(F2_multiplicity_LP(:,j));
    TransEnergy(:,j)                    = Energy';
    
    PhotonTransAtNormalIncDS(:,j)       = PercentTransmissionDS(:,j).*Flux;
    PhotonTransAtNormalIncUS(:,j)       = PercentTransmissionUS(:,j).*Flux;
    
%     PhotonTransAtNormalIncDS(:,j)       = PercentAbsorptionGe(:,j).*PercentTransmissionDS(:,j).*F2_normalized(:,j).*Flux;
%     PhotonTransAtNormalIncUS(:,j)       = PercentAbsorptionGe(:,j).*PercentTransmissionUS(:,j).*F2_normalized(:,j).*Flux;
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
line([0 numhkls], [max(BeamLineFlux(:,1)) max(BeamLineFlux(:,1))], ...
    'LineStyle', '--', 'Color', 'r', 'LineWidth', 2)
title('Diffraction energy - First and the last TOA only')
grid on

subplot(2,3,3)
plot(1:1:numhkls, F2_normalized(:,1), 'k.-')
xlabel('hkl id')
ylabel('F^2*multiplicity*LP (normalized to maximum)')
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
ylabel('transmission through sample (%)')
grid on
axis([1 numhkls 0 100])

subplot(2,3,5)
plot(1:1:numhkls, PercentAbsorptionGe(:,1)*100, 'b.')
hold on
plot(1:1:numhkls, PercentAbsorptionGe(:,end)*100, 'r.')
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('nominal detector efficiency (%)')
grid on
axis([1 numhkls 0 100])

subplot(2,3,6)
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
ylabel('number of photons detected (-)')
grid on

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