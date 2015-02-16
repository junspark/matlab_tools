clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCC Fe
latticeParms_bcc    = 2.87;                        % IN Angstrom
hkls_bcc            = load('bcc.hkls');
d_hkl_bcc           = PlaneSpacings(latticeParms_bcc, 'cubic', hkls_bcc');
pkid_fit            = 4:8;
sqrt_hkls_bcc       = sqrt(sum(hkls_bcc(pkid_fit, :).*hkls_bcc(pkid_fit, :),2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HORZ
ChToEnergyConversion = [0.0926   -0.0754];
TOA  = 6.9678;

lambda_hkl_bcc  = 2.*d_hkl_bcc*sind(TOA/2);

E_hkl_bcc   = Angstrom2keV(lambda_hkl_bcc);

x       = 1:1:2048;
E_grid  = ChToEnergyConversion(1)*x + ChToEnergyConversion(2);

pname_data  = 'C:\Users\parkjs\Box Sync\Projects\Operations\UserSupport\mach_201502_edd\Lap7_1';
dirlist = dir([pname_data, '/horz*']);
data        = zeros(2048,1);
dataStack   = [];
for i = 1:1:length(dirlist)
    pfname  = fullfile(pname_data, dirlist(i).name);
    datai   = load(pfname);
    data    = data + datai;
    
    figure(1)
    plot(E_grid, datai, 'k.')
    hold on
    plot(E_hkl_bcc, ones(length(E_hkl_bcc), 1), 'b^')
    plot(E_hkl_bcc(pkid_fit), ones(length(pkid_fit), 1), 'r^')
    grid on
    xlabel('Energy (keV)')
    ylabel('counts')
    
    for j = 1:1:length(pkid_fit)
        idx1    = find(E_hkl_bcc(pkid_fit(j))-3.5 < E_grid);
        idx2    = find(E_hkl_bcc(pkid_fit(j))+3.5 < E_grid);
        xdata   = E_grid(idx1:idx2);
        ydata   = datai(idx1:idx2);
        
        p0  = [ ...
            max(ydata); ...
            1; ...
            0.5; ...
            E_hkl_bcc(pkid_fit(j)); ...
            0; ...
            0; ...
            ];
        pLB = [ ...
            0; ...
            0; ...
            0; ...
            E_hkl_bcc(pkid_fit(j)) - 3; ...
            -inf; ...
            -inf; ...
            ];
        pUB = [ ...
            inf; ...
            inf; ...
            1; ...
            E_hkl_bcc(pkid_fit(j)) + 3; ...
            inf; ...
            inf; ...
            ];
        
        p   = lsqcurvefit(@pfunc, p0, xdata, ydata, pLB, pUB);
        
        yfit0   = pfunc(p0, xdata);
        yfit    = pfunc(p, xdata);
        
        plot(xdata, ydata, 'r.')
        plot(xdata, yfit0, 'g-')
        plot(xdata, yfit, 'b-')
        Efit(i,j)   = p(4);
    end
    figure(1)
    hold off
    
    lambda      = keV2Angstrom(Efit(i,:));
    d_hkl_fit   = lambda./2/sind(TOA/2);
    
    xdata       = 1./sqrt_hkls_bcc;
    ydata       = d_hkl_fit';
    p           = polyfit(xdata, ydata, 1);
    
    xdata
    ydata
    a_fit(i)    = p(1);
    
    figure(10);
    plot(xdata, ydata, 'o')
    hold on
    plot(xdata, polyval(p, xdata), '-') 
    xlabel('1/sqrt(hh+kk+ll)')
    ylabel('d_{hkl} (Angstrom)')
    pause
end

figure,
plot(E_grid, data, 'k.')
hold on
plot(E_hkl_bcc, ones(length(E_hkl_bcc), 1), 'b^')
plot(E_hkl_bcc(pkid_fit), ones(length(pkid_fit), 1), 'r^')
grid on
xlabel('Energy (keV)')
ylabel('counts')