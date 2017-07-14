clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOA  = 4.9995;
ChToEnergyConversion  = [0.0346714 0.0125024];
MeasurementPlane        = 'v';

% TOA  = 4.90376;
% ChToEnergyConversion  = [0.0341265 -0.00588104];
% MeasurementPlane        = 'h';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCC Fe
latticeParms    = 3.61;                        % IN Angstrom
hkls            = load('fcc.hkls');
d_hkl           = PlaneSpacings(latticeParms, 'cubic', hkls');
pkid_fit        = 1:5;
sqrt_hkls       = sqrt(sum(hkls(pkid_fit, :).*hkls(pkid_fit, :),2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAR FILE DESIGNATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname_pypar     = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/edd_6bm_2017-2/startup_jun17';
fname_pypar     = 'exposure_record_full.par';
pfname_pypar    = fullfile(pname_pypar, fname_pypar);
pardata         = ReadSpecParFile(pfname_pypar, 'Version', 'mpe_standard');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH WHERE DATA FILES LIVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname   = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/edd_6bm_2017-2/startup_jun17';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PEAK FITTING OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkfitting_pars.pfunc_type   = 'pseudovoigt';
pkfitting_pars.pbkg_order   = 2;

determine_a0    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XS  = pardata.samX;
YS  = pardata.samY;
ZS  = pardata.samZ;

numDV   = length(pardata.day);
for i = 1:1:numDV
    fname   = pardata.det1_fname{i};
    
    if ~isempty(strfind(fname, 'afrl_a52_i'))
        pname_data  = fullfile(pname, fname);
        for j = 1:57
            fname_data  = sprintf('%s-%03d-hv.xy', fname, j);
            pfname_data = fullfile(pname, fname, fname_data);
            
            data    = load(pfname_data);
            for k = 1:1:2
                if k == 1
                    MeasurementPlane        = 'v';
                    TOA  = 4.9995;
                    ChToEnergyConversion  = [0.0346714 0.0125024];
                    
                else
                    MeasurementPlane        = 'h';
                    TOA  = 4.90376;
                    ChToEnergyConversion  = [0.0341265 -0.00588104];
                end
                datai= data(:,k);
                
                fname_fit   = sprintf('%s.%s.fit.mat', fname_data, MeasurementPlane);
                pfname_fit  = fullfile(pname_data, fname_fit);
                
                lambda_hkl  = 2.*d_hkl*sind(TOA/2);
                E_hkl       = Angstrom2keV(lambda_hkl);
                
                E_grid      = 1:1:8192;
                E_grid      = ChToEnergyConversion(1)*E_grid + ChToEnergyConversion(2);
                
%                 figure(1)
%                 plot(E_grid, datai, 'k.')
%                 hold on
%                 plot(E_hkl, ones(length(E_hkl), 1), 'b^')
%                 plot(E_hkl(pkid_fit), ones(length(pkid_fit), 1), 'r^')
%                 grid on
%                 xlabel('Energy (keV)')
%                 ylabel('counts')
                
                for jj = 1:1:length(pkid_fit)
                    idx1    = find(E_hkl(pkid_fit(jj))-3.0 < E_grid);
                    idx2    = find(E_hkl(pkid_fit(jj))+3.0 < E_grid);
                    
                    xdata   = E_grid(idx1:idx2);
                    ydata   = datai(idx1:idx2);
                    
                    %%% PV
                    switch lower(pkfitting_pars.pfunc_type)
                        case 'pseudovoigt'
                            p0  = [ max(ydata); 0.5; 0.01;     E_hkl(pkid_fit(jj));     0;   0; ];
                            pLB = [          0;   0;    0; E_hkl(pkid_fit(jj)) - 3; -inf; -inf; ];
                            pUB = [        inf; inf;    1; E_hkl(pkid_fit(jj)) + 3;  inf;  inf; ];
                        case 'splitpseudovoigt'
                            p0  = [ max(ydata); 0.5; 0.5; 0.01; 0.01;     E_hkl(pkid_fit(jj));    0;    0; ];
                            pLB = [          0;   0;   0;    0;    0; E_hkl(pkid_fit(jj)) - 3; -inf; -inf; ];
                            pUB = [        inf; inf; inf;    1;    1; E_hkl(pkid_fit(jj)) + 3; +inf; +inf; ];
                        otherwise
                            disp('Unknown peak function!!')
                            return
                    end
                    pkfitting_pars.xdata    = xdata;
                    [p, rn(jj), ~, ef(jj)]    = lsqcurvefit(@pfunc_switch, p0, pkfitting_pars, ydata, pLB, pUB);
                    
                    yfit0   = pfunc_switch(p0, pkfitting_pars);
                    yfit    = pfunc_switch(p, pkfitting_pars);
                    
%                     plot(xdata, ydata, 'r.')
%                     plot(xdata, yfit0, 'g-')
%                     plot(xdata, yfit, 'b-')
                    
                    pout    = pkfit_MapResult(pkfitting_pars, p);
                    
                    Afit(jj) = pout(1);
                    gfit{jj} = pout(2:3);
                    nfit{jj} = pout(4:5);
                    Efit(jj) = pout(6);
                    bkg{jj}  = p(7:end);
                    
                    Rwp(jj)  = ErrorRwp(ydata, yfit);
                    Re(jj)   = ErrorRe(ydata, yfit);
                    Rp(jj)   = ErrorRp(ydata, yfit);
                end
%                 figure(1)
%                 hold off
                
                if determine_a0 == 1
                    lambda      = keV2Angstrom(Efit);
                    d_hkl_fit   = lambda./2/sind(TOA/2);
                    
                    xdata       = 1./sqrt_hkls;
                    ydata       = d_hkl_fit';
                    
                    [p, a0_fit_s]   = polyfit(xdata, ydata, 1);
                    a0_fit          = p;
                    
                    figure(10);
                    plot(xdata, ydata, 'o')
                    hold on
                    plot(xdata, polyval(p, xdata), '-')
                    xlabel('1/sqrt(hh+kk+ll)')
                    ylabel('d_{hkl} (Angstrom)')
                    title('a0 is the slope if the data are from reference state and material is cubic')
                else
                    a0_fit      = nan;
                    a0_fit_s    = nan;
                end
%                 pause(1)
                save(pfname_fit, 'Afit', 'gfit', 'nfit', 'Efit', 'bkg', 'Rwp', 'Re', 'Rp', 'rn', 'ef', 'a0_fit', 'a0_fit_s');
            end
        end
    end
end
disp(sprintf('fitting completed'))
