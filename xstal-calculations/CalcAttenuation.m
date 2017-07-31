clear all
close all
clc

E_keV   = [0:0.01:130]';
E_MeV   = E_keV./1000; % MeV


[mac, ~, rho]   = PhotonAttenuation('Fe', E_MeV);

x = (1./mac)./rho; % cm

figure(1)
plot(E_keV, x)
xlabel('Energy [keV]')
ylabel('1/(mu*rho) [cm]')

x           = 0.7; % cm
I_by_Io     = exp(-mac*rho*x);
figure(2)
plot(E_keV, I_by_Io)
xlabel('Energy [keV]')
ylabel('I/Io [-]')
