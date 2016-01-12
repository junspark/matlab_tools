clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BeamLineFlux    = load('bm_flux.data');
TakeOffAngle    = 3:0.5:15;                    % IN deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda  = keV2Angstrom(BeamLineFlux(:,1));

for i = 1:1:length(TakeOffAngle)
    d_coverage(:,i) = lambda./(sind(TakeOffAngle(i)/2)*2);
end

figure(1)
semilogy(BeamLineFlux(:,1), d_coverage(:,1), 'b.')
hold on
semilogy(BeamLineFlux(:,1), d_coverage(:,end), 'r.')
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('Energy (keV)')
ylabel('d-spacing (Angstrom)')
grid on

% figure(2)
% semilogy(BeamLineFlux(:,1), d_coverage(:,1), 'b.')
% hold on
% semilogy(BeamLineFlux(:,1), d_coverage(:,end), 'r.')
% legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
% xlabel('Energy (keV)')
% ylabel('d-spacing (Angstrom)')
% grid on