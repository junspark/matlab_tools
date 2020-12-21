function status = FitPeaksPerPolImage(pfname, Material, Analysis_Options, varargin)

% default options
optcell = {...
    'ShowPlot', false, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

%%%%% 
disp('###########################')
fprintf('Fitting peaks in %s\n', pfname)
disp('###########################')

pfname_polimg   = sprintf('%s.polimg.mat', pfname);
pfname_pkfit    = sprintf('%s.pkfit.mat', pfname);
pfname_pkfit_tbl    = sprintf('%s.pkfit.csv', pfname);

polimg  = load(pfname_polimg);

CakePrms    = polimg.cakeprms;
Instr       = polimg.instr;
polimg      = polimg.polimg;

switch opts.ShowPlot
    case {true, 1}
        figure(2)
        subplot(1,2,1)
        imagesc(log(abs(polimg.intensity_in_tth_grid))), axis square tight
        hold off

        subplot(1,2,2)
        plot(polimg.tth_grid, polimg.intensity_in_tth_grid)
        hold off
end

csv_table   = [];
for kkk = 1:1:CakePrms.bins(1)
    csv_table_kkk   = polimg.azimuth(kkk);
    
    fprintf('Looking at azimuthal bin %d of %d\n', kkk, CakePrms.bins(1))
    
    x   = polimg.tth_grid;
    y   = polimg.intensity_in_tth_grid(kkk,:);
    
    switch opts.ShowPlot
        case {true, 1}
            figure(11)
            subplot(1,2,1)
            plot(x, y, 'k.')
            hold on
            plot(tth, mean(y), 'g^')
            axis([min(x) max(x) min(y) max(y)+100])
            xlabel('radial distance (mm)')
            ylabel('intensity (arb. units)')
            title(['bin number : ', num2str(kkk)])
    end
    
    for mmm = 1:1:Material.numbounds
        numpks  = length(Material.pkidx{mmm});
        fprintf('Looking at bound number %d of %d with %d peaks\n', mmm, Material.numbounds, numpks);
        
        pkrange = Material.pkrange(:,mmm);
        
        idx = find(x >= pkrange(1) & x <= pkrange(2));
        xr  = x(idx)';
        yr  = y(idx)';
        
        if kkk == 1
            idx_max = [];
            peakdet_thresh  = 0.5;
            while (length(idx_max) ~= numpks) && (peakdet_thresh > 0)
                [idx_max, idx_min]  = peakdet(yr, peakdet_thresh, xr);
                peakdet_thresh      = peakdet_thresh - 0.1;
            end
            
            %                             pr0 = [];
            %                             if length(idx_max) == numpks
            %                                 for nnn = 1:1:numpks
            %                                 % pr0 =
            %                                 end
            %                             else
            %                                 for nnn = 1:1:numpks
            %                                     % pr0 =
            %                                 end
            %                             end
            
            pkidx   = Material.pkidx{mmm};
            tth     = Material.tth;
            tth_UB  = Material.tth_UB;
            tth_LB  = Material.tth_LB;
            
            pr0     = [];
            pr0_LB  = [];
            pr0_UB  = [];
            for nnn = 1:1:numpks
                pr0 = [ ...
                    pr0;
                    idx_max(nnn,2)/10;
                    0.05;
                    0.5;
                    tth(pkidx(nnn));
                    ];
                
                pr0_LB  = [ ...
                    pr0_LB; ...
                    0; ...
                    0; ...
                    0; ...
                    tth_LB(pkidx(nnn)); ...
                    ];
                
                pr0_UB  = [ ...
                    pr0_UB; ...
                    inf; ...
                    inf; ...
                    1; ...
                    tth_UB(pkidx(nnn)); ...
                    ];
            end
            
            pr0 = [pr0; ...
                0; ...
                yr(1); ...
                ];
            
            pr0_LB  = [pr0_LB; ...
                -inf; ...
                -inf; ...
                ];
            
            pr0_UB  = [pr0_UB; ...
                inf; ...
                inf; ...
                ];
        else
            pr0     = pr_previous_azimuth{mmm};
            pr0_LB  = pr_LB_previous_azimuth{mmm};
            pr0_UB  = pr_UB_previous_azimuth{mmm};
        end
        
        pkpars.pfunc_type   = Analysis_Options.PkFuncOptions.pfunc_type;
        pkpars.pbkg_order   = Analysis_Options.PkFuncOptions.pbkg_order;
        pkpars.xdata        = xr;
        
        [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc_switch, pr0, pkpars, yr, ...
            pr0_LB, pr0_UB, Analysis_Options.PkFitOptimizationOptions);
        y0	= pfunc_switch(pr0, pkpars);
        yf	= pfunc_switch(pr, pkpars);
        
        % [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc, pr0, xr, yr, ...
        %     [], [], Analysis_Options.PkFitOptions);
        % y0  = pfunc(pr0,xr);
        % yf  = pfunc(pr,xr);
        
        %%% UPDATE PREVIOUS AZIMUTH PR FOR NEXT AZIM FITING
        pr_previous_azimuth{mmm}    = pr;
        pr_LB_previous_azimuth{mmm} = pr0_LB;
        pr_UB_previous_azimuth{mmm} = pr0_UB;
        
        switch opts.ShowPlot
            case {true, 1}
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
                xlabel('radial distance (mm)')
                ylabel('intensity (arb. units)')
                title(['bound number : ', num2str(mmm)])
                hold off
        end
        
        %%% MAPPING NEEDS TO BE UPDATED WITH A SWITCH
        % pro  = pkfitResultMapping(pkpars, pr);
        rwp = ErrorRwp(yr, yf);
        
        pkfit.amp(kkk,mmm)  = pr(1);
        pkfit.fwhm(kkk,mmm) = pr(2);
        pkfit.mix(kkk,mmm)  = pr(3);
        pkfit.rho(kkk,mmm)  = pr(4);
        pkfit.bkg{kkk,mmm}  = pr(5:end);
        pkfit.rsn(kkk,mmm)  = rsn;
        pkfit.ef(kkk,mmm)   = ef;
        pkfit.rwp(kkk,mmm)  = rwp;
        
        csv_table_kkk   = [csv_table_kkk pr' rsn ef rwp];
    end
    csv_table   = [ 
        csv_table; ...
        csv_table_kkk;];
    switch opts.ShowPlot
        case {true, 1}
            figure(11)
            subplot(1,2,1)
            hold off
    end    
end

if Analysis_Options.save_fits
    disp('###########################')
    fprintf('Saving peak fits in %s\n', pfname_pkfit)
    save(pfname_pkfit, 'pfname', 'Material', 'pkfit', 'pkpars', 'CakePrms', 'Instr', 'polimg');
    csvwrite(pfname_pkfit_tbl, csv_table);
else
    disp('###########################')
    fprintf('Not saving peak fits for %s\n', pfname)
end
status  = 1;
