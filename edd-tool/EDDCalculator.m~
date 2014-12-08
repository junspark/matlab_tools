clear all
close all
clc

MaterialName    = 'Al';
latticeParms    = 4.050;                % IN Angstrom
SampleThickness = 0.25;                 % IN cm
hkls            = load('fcc.hkls');
d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

BeamLineFlux    = load('bm_flux.data');
TakeOffAngle    = 7:1:13;                    % IN deg

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
subplot(2,2,1)
plot(TransEnergy(:,1), '.')
hold on
plot(TransEnergy(:,end), 'r.')
legend('7', '13')
xlabel('hkl id')
ylabel('Energy (keV)')
grid on

subplot(2,2,2)
semilogy(PhotonTransAtNormalIncidence(:,1), 'b.')
hold on
semilogy(PhotonTransAtNormalIncidence(:,end), 'r.')
legend('7', '13')
xlabel('hkl id')
ylabel('Number of photons transmitted')
grid on

subplot(2,2,3)
semilogy(BeamLineFlux(:,1), BeamLineFlux(:,2), 'k.')
xlabel('Energy (keV)')
ylabel('Number of photons')
grid on

subplot(2,2,4)
plot(PercentTransmission(:,1), 'b.')
hold on
plot(PercentTransmission(:,end), 'r.')
legend('7', '13')
xlabel('hkl id')
ylabel('Percent transmission')
grid on