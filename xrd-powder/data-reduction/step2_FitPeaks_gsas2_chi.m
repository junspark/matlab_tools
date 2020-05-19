clear all
close all
clc

addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%% DATA REDUCTION FLAGS
Analysis_Options.PkFitOptimizationOptions   = optimset(...
    'MaxIter', 5e5, ...
    'MaxFunEvals',3e5);

Analysis_Options.make_polimg    = 0;
Analysis_Options.save_polimg    = 0;
Analysis_Options.fits_spectra   = 1;
Analysis_Options.save_fits      = 1;
Analysis_Options.find_instrpars = 0;
Analysis_Options.save_instrpars = 0;
Analysis_Options.find_detpars	= 0;
Analysis_Options.generateESG    = 0;

%%% INPUT PARAMETERS
XRDIMAGE.ExpID              = 'mli_dec19';
XRDIMAGE.MetaDataFile       = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\mli_nov19_Tomo.par';
XRDIMAGE.Image.pname_chi    = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\data\dexela\corrected';
XRDIMAGE.Image.pname        = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\data\dexela\corrected';
XRDIMAGE.Image.fbase        = 'am316_ss_id105_1';
XRDIMAGE.Image.fnumber      = 324:326;
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.corrected    = 1;   % 0 - UNCORRECTED, 1 = SUM, 2 = COR32
XRDIMAGE.Image.dettype      = 6;   % 3 = GE3, 5 = GE5, 1234 = HYDRA1234, 6 = DEX IN C

%%% FITTING PARAMETER
XRDIMAGE.pkListFile         = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\waxs\am316_ss_id105_1.pkslst';

%%% ANALYSIS STARTS HERE
metadata    = ReadSpecParFile(XRDIMAGE.MetaDataFile, 'Version', 'mli_nov19_c');

%%% READ PEAK LIST FROM GSAS2
fid     = fopen(XRDIMAGE.pkListFile, 'r');
pklist  = [];
while ~feof(fid)
    linedata    = fgetl(fid);
    switch linedata(1)
        case '['
            pklist  = [pklist; eval(linedata)];
    end
end
fclose(fid);
numpks  = size(pklist ,1);
switch XRDIMAGE.Image.dettype
    case 6
        disp(sprintf('%d = DEXELA', XRDIMAGE.Image.dettype));
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            
            idx_fbase       = ismember(metadata.det1_fname, XRDIMAGE.Image.fbase);
            idx_fnumber     = ismember(metadata.det1_fnum, XRDIMAGE.Image.fnumber(iii));
            idx_metadata    = find(idx_fbase & idx_fnumber);
            
            switch ~isempty(idx_metadata) && (length(idx_metadata) == 1)
                case true
                    fext   = 'tif';
                    
                    fimctrl = sprintf('%s_%s_%s.imctrl', XRDIMAGE.ExpID, fext, XRDIMAGE.Image.fbase);
                    imctrl  = ReadGSAS2Imctrl(fimctrl);
                    
                    dazm        = (imctrl.azm_range(2) - imctrl.azm_range(1))/imctrl.azm_steps;
                    azm_grid    = linspace(imctrl.azm_range(1)+dazm/2, imctrl.azm_range(2)-dazm/2, imctrl.azm_steps);
                    
                    idx_azm             = azm_grid >= 360;
                    azm_grid(idx_azm)   = azm_grid(idx_azm) - 360;
                    switch XRDIMAGE.Image.corrected
                        case 1
                            fname_pattern   = sprintf('%s_bkg_corrected_%s', froot, fext);
                            fname_pkfit     = sprintf('%s.pkfit.mat',fname_pattern);
                            pfname_pkfit    = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fname_pkfit);
                            
                            %%% INIT FOR PARFOR
                            pkfit_amp   = zeros(imctrl.azm_steps, numpks);
                            pkfit_mix	= zeros(imctrl.azm_steps, numpks);
                            pkfit_fwhm	= zeros(imctrl.azm_steps, numpks);
                            pkfit_rho	= zeros(imctrl.azm_steps, numpks);
                            pkfit_bkg	= cell(imctrl.azm_steps, numpks);
                            pkfit_rsn	= zeros(imctrl.azm_steps, numpks);
                            pkfit_ef	= zeros(imctrl.azm_steps, numpks);
                            pkfit_rwp	= zeros(imctrl.azm_steps, numpks);
                            pkfit_xr    = cell(imctrl.azm_steps, numpks);
                            pkfit_yr    = cell(imctrl.azm_steps, numpks);
                            pkfit_eta   = azm_grid';
                            pkfit_samX  = metadata.samX(idx_metadata);
                            pkfit_samY  = metadata.samY(idx_metadata);
                            pkfit_samZ  = metadata.samZ(idx_metadata);
                            pkfit_rX    = metadata.aX(idx_metadata);
                            pkfit_rY    = metadata.aY(idx_metadata);
                            pkfit_rZ    = metadata.aZ(idx_metadata);
                            %                             pkfit_samX2 = metadata.samX2(idx_metadata);
                            %                             pkfit_samY2 = metadata.samY2(idx_metadata);
                            %                             pkfit_samZ2 = metadata.samZ2(idx_metadata);
                            %                             pkfit_samP  = metadata.samOther(idx_metadata);
                            %                             pkfit_enc1  = metadata.encoder1(idx_metadata);
                            %                             pkfit_enc2  = metadata.encoder2(idx_metadata);
                            %                             pkfit_enc3  = metadata.encoder3(idx_metadata);
                            %                             pkfit_enc4  = metadata.encoder4(idx_metadata);
                            %                             pkfit_enc5  = metadata.encoder5(idx_metadata);
                            %                             pkfit_enc6  = metadata.encoder6(idx_metadata);
                            %                             pkfit_enc7  = metadata.encoder7(idx_metadata);
                            %                             pkfit_enc8  = metadata.encoder8(idx_metadata);
                            %                             pkfit_enc9  = metadata.encoder9(idx_metadata);
                            %                             pkfit_enc10 = metadata.encoder10(idx_metadata);
                            %                             pkfit_ev1   = metadata.ev1(idx_metadata);
                            %                             pkfit_ev2   = metadata.ev2(idx_metadata);
                            %                             pkfit_ev3   = metadata.ev3(idx_metadata);
                            %                             pkfit_ev4   = metadata.ev4(idx_metadata);
                            %                             pkfit_ev5   = metadata.ev5(idx_metadata);
                            %                             pkfit_ev6   = metadata.ev6(idx_metadata);
                            %                             pkfit_ev7   = metadata.ev7(idx_metadata);
                            %                             pkfit_ev8   = metadata.ev8(idx_metadata);
                            %                             pkfit_ev9   = metadata.ev9(idx_metadata);
                            %                             pkfit_ev10  = metadata.ev10(idx_metadata);
                            
                            for kkk = 1:imctrl.azm_steps
                                fchi    = sprintf('%s_A%d.chi', fname_pattern, fix(azm_grid(kkk)));
                                pfchi   = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fchi);
                                switch isfile(pfchi)
                                    case true
                                        disp(sprintf('%s DOES EXISTS\n', fchi));
                                        I_vs_tth    = table2array(readtable(pfchi, 'FileType', 'text', ...
                                            'ReadVariableNames', false, ...
                                            'HeaderLines', 4));
                                        
                                        x   = I_vs_tth(:,1);
                                        y   = I_vs_tth(:,2);
                                        
                                        figure(11)
                                        subplot(1,2,1)
                                        plot(x, y, 'k.')
                                        hold on
                                        plot(pklist(:,1), mean(y), 'g^')
                                        axis([min(x) max(x) 0 max(y)+100])
                                        xlabel('2theta (\theta)')
                                        ylabel('intensity (arb. units)')
                                        title(['bin number : ', num2str(kkk)])
                                        for mmm = 1:1:numpks
                                            xri = pklist(mmm, 1)-0.15;
                                            xrf = pklist(mmm, 1)+0.15;
                                            idx_metadata = find(x >= xri & x <= xrf);
                                            
                                            xr  = x(idx_metadata);
                                            yr  = y(idx_metadata);
                                            
                                            pkfit_xr{kkk,mmm}   = xr;
                                            pkfit_yr{kkk,mmm}   = yr;
                                            
                                            yr_nan  = sum(yr == 0);
                                            switch yr_nan > 1
                                                case true
                                                    warning('there are NAN in the peak profile; will not fit.')
                                                    pkfit_amp(kkk,mmm)  = nan;
                                                    pkfit_mix(kkk,mmm)  = nan;
                                                    pkfit_fwhm(kkk,mmm) = nan;
                                                    pkfit_rho(kkk,mmm)  = nan;
                                                    pkfit_bkg{kkk,mmm}  = [nan nan];
                                                    pkfit_rsn(kkk,mmm)  = nan;
                                                    pkfit_ef(kkk,mmm)   = nan;
                                                    pkfit_rwp(kkk,mmm)  = nan;
                                                otherwise
                                                    disp('fitting peak.')
                                                    pr0 = pklist(mmm, 1:2:end)';
                                                    pr0 = [...
                                                        max(yr)/5 ...
                                                        0.1 ...
                                                        0.5 ...
                                                        pr0(1) ...
                                                        0 ...
                                                        yr(1)];
                                                    prLB    = [0 0 0 xri -inf -inf]';
                                                    prUB    = [inf inf 1 xrf inf inf]';
                                                    
                                                    [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc, pr0, xr, yr, ...
                                                        prLB, prUB, Analysis_Options.PkFitOptimizationOptions);
                                                    
                                                    y0	= pfunc(pr0, xr);
                                                    yf	= pfunc(pr, xr);
                                                    
                                                    figure(11)
                                                    subplot(1,2,1)
                                                    plot(xr, yr, 'b.')
                                                    plot(xr, y0, 'r-')
                                                    plot(xr, yf, 'g-')
                                                    
                                                    subplot(1,2,2)
                                                    plot(xr, yr, 'b.')
                                                    hold on
                                                    plot(xr, y0, 'r-')
                                                    plot(xr, yf, 'g-')
                                                    xlabel('2theta (\theta)')
                                                    ylabel('intensity (arb. units)')
                                                    title(['peak number : ', num2str(mmm)])
                                                    hold off
                                                    
                                                    pkfit_amp(kkk,mmm)  = pr(1);
                                                    pkfit_mix(kkk,mmm)  = pr(2);
                                                    pkfit_fwhm(kkk,mmm) = pr(3);
                                                    pkfit_rho(kkk,mmm)  = pr(4);
                                                    pkfit_bkg{kkk,mmm}  = pr(5:end);
                                                    pkfit_rsn(kkk,mmm)  = rsn;
                                                    pkfit_ef(kkk,mmm)   = ef;
                                                    pkfit_rwp(kkk,mmm)  = ErrorRwp(yr, yf);
                                            end
                                        end
                                    otherwise
                                        error('%s DOES NOT EXIST\n', pfchi);
                                end
                            end
                            
                            switch Analysis_Options.save_fits
                                case true
                                    disp('###########################')
                                    fprintf('Saving peak fits in %s\n', pfname_pkfit)
                                    save(pfname_pkfit, 'pkfit*')
                                otherwise
                                    disp('###########################')
                                    fprintf('Not saving peak fits for %s\n', pfname_pkfit)
                            end
                        otherwise
                            warning('XRDIMAGE.Image.corrected = %d is not a valid option', XRDIMAGE.Image.corrected);
                            return
                    end
                otherwise
                    warning('%s metadata needs investigation')
            end
        end
    case 1234
        disp(sprintf('%d = HYDRA1234', XRDIMAGE.Image.dettype));
        
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            
            idx_fbase       = ismember(metadata.det1_fname, XRDIMAGE.Image.fbase);
            idx_fnumber     = ismember(metadata.det1_fnum, XRDIMAGE.Image.fnumber(iii));
            idx_metadata    = find(idx_fbase & idx_fnumber);
            switch ~isempty(idx_metadata) && (length(idx_metadata) == 1)
                case true
                    for jjj = 1:1:4
                        fext    = sprintf('ge%d', jjj);
                        fimctrl = sprintf('%s_%s_%s.imctrl', XRDIMAGE.ExpID, fext, XRDIMAGE.Image.fbase);
                        imctrl  = readtable(fimctrl, 'FileType', 'text', 'ReadVariableNames', false, ...
                            'delimiter', ':');
                        azm_range   = eval(cell2mat(imctrl.Var2(find(ismember(imctrl.Var1, 'LRazimuth')))));
                        azm_offset  = str2double(cell2mat(imctrl.Var2(find(ismember(imctrl.Var1, 'azmthOff')))));
                        azm_steps   = str2double(cell2mat(imctrl.Var2(find(ismember(imctrl.Var1, 'outAzimuths')))));
                        tth_range   = eval(cell2mat(imctrl.Var2(find(ismember(imctrl.Var1, 'IOtth')))));
                        tth_steps   = str2double(cell2mat(imctrl.Var2(find(ismember(imctrl.Var1, 'outChannels')))));
                        
                        dazm        = (azm_range(2) - azm_range(1))/azm_steps;
                        azm_grid    = linspace(azm_range(1)+dazm/2, azm_range(2)-dazm/2, azm_steps);
                        
                        idx_azm     = azm_grid >= 360;
                        azm_grid(idx_azm)   = azm_grid(idx_azm) - 360;
                        
                        switch XRDIMAGE.Image.corrected
                            case 1
                                fname_pattern   = sprintf('%s.%s.sum', froot, fext);
                                
                                fname_pkfit    = sprintf('%s.pkfit.mat',fname_pattern);
                                pfname_pkfit    = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fname_pkfit);
                                
                                %%% INIT FOR PARFOR
                                pkfit_amp   = zeros(azm_steps, numpks);
                                pkfit_mix	= zeros(azm_steps, numpks);
                                pkfit_fwhm	= zeros(azm_steps, numpks);
                                pkfit_rho	= zeros(azm_steps, numpks);
                                pkfit_bkg	= cell(azm_steps, numpks);
                                pkfit_rsn	= zeros(azm_steps, numpks);
                                pkfit_ef	= zeros(azm_steps, numpks);
                                pkfit_rwp	= zeros(azm_steps, numpks);
                                pkfit_xr    = cell(azm_steps, numpks);
                                pkfit_yr    = cell(azm_steps, numpks);
                                pkfit_eta   = azm_grid';
                                pkfit_samX  = metadata.samX(idx_metadata);
                                pkfit_samY  = metadata.samY(idx_metadata);
                                pkfit_samZ  = metadata.samZ(idx_metadata);
                                pkfit_rX    = metadata.aX(idx_metadata);
                                pkfit_rY    = metadata.aY(idx_metadata);
                                pkfit_rZ    = metadata.aZ(idx_metadata);
                                pkfit_samX2 = metadata.samX2(idx_metadata);
                                pkfit_samY2 = metadata.samY2(idx_metadata);
                                pkfit_samZ2 = metadata.samZ2(idx_metadata);
                                pkfit_samP  = metadata.samOther(idx_metadata);
                                pkfit_enc1  = metadata.encoder1(idx_metadata);
                                pkfit_enc2  = metadata.encoder2(idx_metadata);
                                pkfit_enc3  = metadata.encoder3(idx_metadata);
                                pkfit_enc4  = metadata.encoder4(idx_metadata);
                                pkfit_enc5  = metadata.encoder5(idx_metadata);
                                pkfit_enc6  = metadata.encoder6(idx_metadata);
                                pkfit_enc7  = metadata.encoder7(idx_metadata);
                                pkfit_enc8  = metadata.encoder8(idx_metadata);
                                pkfit_enc9  = metadata.encoder9(idx_metadata);
                                pkfit_enc10 = metadata.encoder10(idx_metadata);
                                pkfit_ev1   = metadata.ev1(idx_metadata);
                                pkfit_ev2   = metadata.ev2(idx_metadata);
                                pkfit_ev3   = metadata.ev3(idx_metadata);
                                pkfit_ev4   = metadata.ev4(idx_metadata);
                                pkfit_ev5   = metadata.ev5(idx_metadata);
                                pkfit_ev6   = metadata.ev6(idx_metadata);
                                pkfit_ev7   = metadata.ev7(idx_metadata);
                                pkfit_ev8   = metadata.ev8(idx_metadata);
                                pkfit_ev9   = metadata.ev9(idx_metadata);
                                pkfit_ev10  = metadata.ev10(idx_metadata);
                                
                                %%% THIS IS WHERE PARFOR CAN COME
                                parfor kkk = 1:azm_steps
                                    fchi    = sprintf('%s_Azm=_%d.chi', fname_pattern, fix(azm_grid(kkk)));
                                    pfchi   = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fchi);
                                    
                                    switch isfile(pfchi)
                                        case true
                                            disp(sprintf('%s DOES EXISTS\n', fchi));
                                            I_vs_tth    = table2array(readtable(pfchi, 'FileType', 'text', ...
                                                'ReadVariableNames', false, ...
                                                'HeaderLines', 4));
                                            
                                            x   = I_vs_tth(:,1);
                                            y   = I_vs_tth(:,2);
                                            
                                            %                                     figure(11)
                                            %                                     subplot(1,2,1)
                                            %                                     plot(x, y, 'k.')
                                            %                                     hold on
                                            %                                     plot(pklist(:,1), mean(y), 'g^')
                                            %                                     axis([min(x) max(x) 0 max(y)+100])
                                            %                                     xlabel('2theta (\theta)')
                                            %                                     ylabel('intensity (arb. units)')
                                            %                                     title(['bin number : ', num2str(kkk)])
                                            
                                            for mmm = 1:1:numpks
                                                xri = pklist(mmm, 1)-0.15;
                                                xrf = pklist(mmm, 1)+0.15;
                                                idx_metadata = find(x >= xri & x <= xrf);
                                                
                                                xr  = x(idx_metadata);
                                                yr  = y(idx_metadata);
                                                
                                                pkfit_xr{kkk,mmm}   = xr;
                                                pkfit_yr{kkk,mmm}   = yr;
                                                
                                                yr_nan  = sum(yr == 0);
                                                switch yr_nan > 1
                                                    case true
                                                        warning('there are NAN in the peak profile; will not fit.')
                                                        pkfit_amp(kkk,mmm)  = nan;
                                                        pkfit_mix(kkk,mmm)  = nan;
                                                        pkfit_fwhm(kkk,mmm) = nan;
                                                        pkfit_rho(kkk,mmm)  = nan;
                                                        pkfit_bkg{kkk,mmm}  = [nan nan];
                                                        pkfit_rsn(kkk,mmm)  = nan;
                                                        pkfit_ef(kkk,mmm)   = nan;
                                                        pkfit_rwp(kkk,mmm)  = nan;
                                                    otherwise
                                                        disp('fitting peak.')
                                                        pr0 = pklist(mmm, 1:2:end)';
                                                        pr0 = [...
                                                            max(yr)/5 ...
                                                            0.1 ...
                                                            0.5 ...
                                                            pr0(1) ...
                                                            0 ...
                                                            yr(1)];
                                                        prLB    = [0 0 0 xri -inf -inf]';
                                                        prUB    = [inf inf 1 xrf inf inf]';
                                                        
                                                        [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc, pr0, xr, yr, ...
                                                            prLB, prUB, Analysis_Options.PkFitOptimizationOptions);
                                                        
                                                        y0	= pfunc(pr0, xr);
                                                        yf	= pfunc(pr, xr);
                                                        
                                                        %                                                 figure(11)
                                                        %                                                 subplot(1,2,1)
                                                        %                                                 plot(xr, yr, 'b.')
                                                        %                                                 plot(xr, y0, 'r-')
                                                        %                                                 plot(xr, yf, 'g-')
                                                        %
                                                        %                                                 subplot(1,2,2)
                                                        %                                                 plot(xr, yr, 'b.')
                                                        %                                                 hold on
                                                        %                                                 plot(xr, y0, 'r-')
                                                        %                                                 plot(xr, yf, 'g-')
                                                        %                                                 xlabel('2theta (\theta)')
                                                        %                                                 ylabel('intensity (arb. units)')
                                                        %                                                 title(['peak number : ', num2str(mmm)])
                                                        %                                                 hold off
                                                        
                                                        pkfit_amp(kkk,mmm)  = pr(1);
                                                        pkfit_mix(kkk,mmm)  = pr(2);
                                                        pkfit_fwhm(kkk,mmm) = pr(3);
                                                        pkfit_rho(kkk,mmm)  = pr(4);
                                                        pkfit_bkg{kkk,mmm}  = pr(5:end);
                                                        pkfit_rsn(kkk,mmm)  = rsn;
                                                        pkfit_ef(kkk,mmm)   = ef;
                                                        pkfit_rwp(kkk,mmm)  = ErrorRwp(yr, yf);
                                                end
                                            end
                                            
                                        otherwise
                                            error('%s DOES NOT EXIST\n', pfchi);
                                    end
                                end
                                
                                switch Analysis_Options.save_fits
                                    case true
                                        disp('###########################')
                                        fprintf('Saving peak fits in %s\n', pfname_pkfit)
                                        save(pfname_pkfit, 'pkfit*')
                                    otherwise
                                        disp('###########################')
                                        fprintf('Not saving peak fits for %s\n', pfname_pkfit)
                                end
                            otherwise
                                warning('XRDIMAGE.Image.corrected = %d is not a valid option', XRDIMAGE.Image.corrected);
                                return
                        end
                    end
                otherwise
                    warning('%s metadata needs investigation')
            end
        end
    otherwise
        warning('%d is not a valid option', XRDIMAGE.Image.dettype);
end
return