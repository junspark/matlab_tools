clear all
close all
clc

%%% PEAK TO FIT
d_pk    = [3.53 2.73];
d_UB    = [3.65 3.04];
d_LB    = [3.21 2.54];

edf_cal = Read_s6bmb_edf('./ALU_031320_0001.edf');
for i = 1:1:10
    b(i)    = edf_cal(i).polynomial(1);
    m(i)    = edf_cal(i).polynomial(2);
    tth(i)  = edf_cal(i).tth(1);
end
proot   = './data';
froot   = 'd36s_longscan000';

pfroot  = fullfile(proot, froot, sprintf('%s.xy', froot));
edd     = readedd5_6bm(pfroot);

numscans    = length(edd);
for iii = 1:1:numscans % 1:1:numscans
    edd_iii = edd(iii);
    
    numsteps    = edd_iii.motor_start_end_numstep(3);
    numdets     = length(edd_iii.data);
    
    samXE   = zeros(numsteps,1);
    samYE   = zeros(numsteps,1);
    samZE   = zeros(numsteps,1);
    for jjj = 1:1:numsteps
        samXE(jjj)  = edd_iii.motorpos_all(jjj,1);
        samYE(jjj)  = edd_iii.motorpos_all(jjj,2);
        samZE(jjj)  = edd_iii.motorpos_all(jjj,3);
        
        for kkk = 1:1:numdets
            ydata   = edd_iii.data{kkk}(jjj,:);
            xdata   = 1:1:length(ydata);
            xgrid   = m(kkk)*xdata + b(kkk);
            
            figure(1)
            plot(xgrid, ydata)
            hold on
            
            for mmm = 1:1:length(d_pk)
                Epk     = Angstrom2keV(2*d_pk(mmm)*sind(tth(kkk)/2));
                ELB     = Angstrom2keV(2*d_UB(mmm)*sind(tth(kkk)/2));
                EUB     = Angstrom2keV(2*d_LB(mmm)*sind(tth(kkk)/2));
                
                idx     = (xgrid > ELB) & (xgrid < EUB);
                
                xfit_data   = xgrid(idx)';
                yfit_data   = double(ydata(idx))';
                
                pr0 = [ ...
                    max(yfit_data)/5 ...
                    0.2 ...
                    0.5 ...
                    Epk ...
                    ];
                
                prUB    = [ ...
                    inf ...
                    inf ...
                    1 ...
                    EUB ...
                    ];
                
                prLB    = [ ...
                    0 ...
                    0 ...
                    0 ...
                    ELB ...
                    ];
                
                pr0 = [pr0 ...
                    0 ...
                    50 ...
                    ];
                
                prUB = [prUB ...
                    inf ...
                    inf ...
                    ];
                
                prLB = [prLB ...
                    -inf ...
                    -inf ...
                    ];
                
                pr0     = pr0';
                prUB    = prUB';
                prLB    = prLB';
                
                [pr, rsn{kkk,mmm}(jjj), ~, ef{kkk,mmm}(jjj)]    = lsqcurvefit(@pfunc, pr0, xfit_data, yfit_data, ...
                    prLB, prUB);
                
                yfit0   = pfunc(pr0, xfit_data);
                yfit    = pfunc(pr, xfit_data);
                
                plot(Epk, 20, 'k^')
                plot(ELB, 20, 'r^')
                plot(EUB, 20, 'rv')
                plot(xfit_data, yfit_data, 'r.')
                plot(xfit_data, yfit0, 'g-')
                plot(xfit_data, yfit, 'k-')
                
                PkInt{kkk, mmm}(jjj)    = trapz(xfit_data, yfit);
                BkgInt{kkk, mmm}(jjj)   = diff(polyval(polyint(pr(5:end)'), [xfit_data(1), xfit_data(end)]));
                
                Aout{kkk, mmm}(jjj) = pr(1);
                Gout{kkk, mmm}(jjj) = pr(2);
                nout{kkk, mmm}(jjj) = pr(3);
                Eout{kkk, mmm}(jjj) = pr(4);
                bkgout{kkk, mmm}(jjj,:) = pr(5:end);
            end
            
            hold off
%             return
            pause(0.25)
        end
    end
    save(sprintf('%s_scan%d_fits.mat', froot, iii), ...
        'samXE', 'samYE', 'samZE', ...
        'Aout', 'Gout', 'nout', 'Eout', 'bkgout', 'rsn', 'ef', 'PkInt', 'BkgInt')
end