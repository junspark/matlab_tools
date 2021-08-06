function status = parallel_FitPeaksPerPolImage(fname_pattern, pname_polimg, pname_pkfit, ...
    Material, Analysis_Options, varargin)

% default options
optcell = {...
    'ShowPlot', false, ...
    'MetaDataFieldName', nan, ...
    'MetaDataFieldValues', nan, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

%%%% LOOK AT THIS LATER
if length(opts.MetaDataFieldName) ~= length(opts.MetaDataFieldValues)
    error('metadata field name length and field value length are different')
    status = -1;
elseif ~isnan(opts.MetaDataFieldName) && ~isnan(opts.MetaDataFieldValues)
    for iii = 1:1:length(opts.MetaDataFieldName)
        csv_header{iii} = opts.MetaDataFieldName{iii};
    end
end

%%%%% 
disp('###########################')
fprintf('Fitting peaks in %s\n', pname_polimg)
disp('###########################')

%%% POLIMG FILE TO LOOK AT
fname_polimg    = sprintf('%s.polimg.mat', fname_pattern);
pfname_polimg   = fullfile(pname_polimg, fname_polimg);

%%% PKFIT FILE NAME TO SAVE
fname_pkfit     = sprintf('%s.pkfit.mat', fname_pattern);
pfname_pkfit    = fullfile(pname_pkfit, fname_pkfit);

%%% PKFIT TABLE FILE NAME TO SAVE
fname_pkfit_tbl     = sprintf('%s.pkfit.csv', fname_pattern);
pfname_pkfit_tbl    = fullfile(pname_pkfit, fname_pkfit_tbl);

%%% LOAD POLIMG
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

ct_hdr_items    = 1;
if ~isnan(opts.MetaDataFieldValues)
    for iii = 1:1:length(opts.MetaDataFieldName)
        csv_table_hdr{ct_hdr_items} = opts.MetaDataFieldName{iii};
        ct_hdr_items    = ct_hdr_items + 1;
    end
end

csv_table       = [];
for kkk = 1:1:CakePrms.bins(1)
    fprintf('Looking at azimuthal bin %d of %d\n', kkk, CakePrms.bins(1));
    csv_table_kkk       = polimg.azimuth(kkk);
    if kkk == 1
        csv_table_hdr{ct_hdr_items} = 'eta';
        ct_hdr_items                = ct_hdr_items + 1;
    end
    
    x   = polimg.tth_grid;
    y   = polimg.intensity_in_tth_grid(kkk,:);
    
    switch opts.ShowPlot
        case {true, 1}
            figure(11)
            subplot(1,2,1)
            plot(x, y, 'k.')
            hold on
            axis([min(x) max(x) min(y) max(y)+100])
            xlabel('radial distance (mm)')
            ylabel('intensity (arb. units)')
            title(['bin number : ', num2str(kkk)])
    end
    
    ct_pk   = 0;
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
                [idx_max, ~]    = peakdet(yr, peakdet_thresh, xr);
                peakdet_thresh  = peakdet_thresh - 0.1;
            end
            
            pkidx   = Material.pkidx{mmm};
            tth     = Material.tth;
            tth_UB  = Material.tth_UB;
            tth_LB  = Material.tth_LB;
            
            pr0     = [];
            pr0_LB  = [];
            pr0_UB  = [];
            
            for nnn = 1:1:numpks
                if numpks == 1
                    pr0 = [ ...
                        pr0;
                        max(yr)/20;
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
                else
                    pr0 = [ ...
                        pr0;
                        idx_max(nnn,2)/2;
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
        
%         [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc_switch, pr0, pkpars, yr, ...
%             pr0_LB, pr0_UB, Analysis_Options.PkFitOptimizationOptions);
%         y0	= pfunc_switch(pr0, pkpars);
%         yf	= pfunc_switch(pr, pkpars);
        
        [pr, rsn, resd, ef, ~, ~, jacobian] = lsqcurvefit(@pfunc, pr0, xr, yr, ...
            pr0_LB, pr0_UB, Analysis_Options.PkFitOptimizationOptions);
        y0  = pfunc(pr0,xr);
        yf  = pfunc(pr,xr);
        
        [conf, var]  = confint(pr, resd, jacobian); % variance of fitted params
        var = full(var);
        
        
        %%% UPDATE PREVIOUS AZIMUTH PR FOR NEXT AZIM FITING
        pr_previous_azimuth{mmm}    = pr;
        pr_LB_previous_azimuth{mmm} = pr0_LB;
        pr_UB_previous_azimuth{mmm} = pr0_UB;
        
        switch opts.ShowPlot
            case {true, 1}
                figure(11)
                subplot(1,2,1)
                plot(tth, mean(y), 'g^')
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
                pause
        end
        
        %%% MAPPING NEEDS TO BE UPDATED WITH A SWITCH
        % pro  = pkfitResultMapping(pkpars, pr);
        rwp = ErrorRwp(yr, yf);
        
        pr_pk   = pr(1:end-2);
        pr_bkg  = pr(end-1:end);
        for nnn = 1:1:numpks
            idx_pk_LB   = 1 + (nnn-1)*4;
            idx_pk_UB   = nnn*4;
            
            % integral(@pfunc_switch(pr_one_pk, pkpars), xr(1), xr(end))
            pr_one_pk   = [pr_pk(idx_pk_LB:idx_pk_UB); pr_bkg];
            yf_one_pk   = pfunc_switch(pr_one_pk, pkpars);
            
            q_bkg       = polyint(pr_bkg');
            
            integrated_intensity_with_background    = trapz(xr, yf_one_pk);
            integrated_background                   = diff(polyval(q_bkg,[xr(1) xr(end)]));
            
            integrated_intensity(nnn)   = integrated_intensity_with_background - integrated_background;
            
            % pr_pk(idx_pk_LB:idx_pk_UB)
            % integrated_intensity(nnn)
            csv_table_kkk   = [csv_table_kkk pr_pk(idx_pk_LB:idx_pk_UB)' integrated_intensity(nnn)];
        end
        
        pkfit.Afit{kkk,mmm}                 = pr_pk(1:4:end);
        pkfit.Afit_eb{kkk,mmm}              = sqrt(var(1:4:end));
        pkfit.Afit_95conf_range{kkk,mmm}    = conf(1:4:end,2) - conf(1:4:end,1);
        
        pkfit.gfit{kkk,mmm}                 = pr_pk(2:4:end);
        pkfit.gfit_eb{kkk,mmm}              = sqrt(var(2:4:end));
        pkfit.gfit_95conf_range{kkk,mmm}    = conf(2:4:end,2) - conf(2:4:end,1);
        
        pkfit.nfit{kkk,mmm}                 = pr_pk(3:4:end);
        pkfit.nfit_eb{kkk,mmm}              = sqrt(var(3:4:end));
        pkfit.nfit_95conf_range{kkk,mmm}    = conf(3:4:end,2) - conf(3:4:end,1);
        
        pkfit.cenfit{kkk,mmm}               = pr_pk(4:4:end);
        pkfit.cenfit_eb{kkk,mmm}            = sqrt(var(4:4:end));
        pkfit.cenfit_95conf_range{kkk,mmm}  = conf(4:4:end,2) - conf(4:4:end,1);
        
        pkfit.bkg{kkk,mmm}                  = pr_bkg;
        pkfit.bkg_eb{kkk,mmm}               = sqrt(var(end-1:end));
        pkfit.bkg_95conf_range{kkk,mmm}     = conf((length(pr_pk)+1:length(pr_pk)+2), 2) - ...
            conf((length(pr_pk)+1:length(pr_pk)+2), 1);
        
        pkfit.rsn(kkk,mmm)      = rsn;
        pkfit.ef(kkk,mmm)       = ef;
        pkfit.rwp(kkk,mmm)      = rwp;
        pkfit.I{kkk,mmm}        = integrated_intensity;
        
        csv_table_kkk   = [csv_table_kkk pr_bkg' rsn ef rwp];
        
        if kkk == 1
            for nnn = 1:1:numpks
                ct_pk   = ct_pk + 1;
                
                csv_table_hdr{ct_hdr_items} = sprintf('pk%d-amp', ct_pk);
                ct_hdr_items    = ct_hdr_items + 1;
                
                csv_table_hdr{ct_hdr_items} = sprintf('pk%d-fwhm', ct_pk);
                ct_hdr_items    = ct_hdr_items + 1;
                
                csv_table_hdr{ct_hdr_items} = sprintf('pk%d-mix', ct_pk);
                ct_hdr_items    = ct_hdr_items + 1;
                
                csv_table_hdr{ct_hdr_items} = sprintf('pk%d-pos', ct_pk);
                ct_hdr_items    = ct_hdr_items + 1;
                
                csv_table_hdr{ct_hdr_items} = sprintf('pk%d-I', ct_pk);
                ct_hdr_items    = ct_hdr_items + 1;
            end
            
            csv_table_hdr{ct_hdr_items} = sprintf('bnd%d-bkg0', mmm);
            ct_hdr_items    = ct_hdr_items + 1;
            
            csv_table_hdr{ct_hdr_items} = sprintf('bnd%d-bkg1', mmm);
            ct_hdr_items    = ct_hdr_items + 1;
            
            csv_table_hdr{ct_hdr_items} = sprintf('bnd%d-rsn', mmm);
            ct_hdr_items    = ct_hdr_items + 1;
            
            csv_table_hdr{ct_hdr_items} = sprintf('bnd%d-ef', mmm);
            ct_hdr_items    = ct_hdr_items + 1;
            
            csv_table_hdr{ct_hdr_items} = sprintf('bnd%d-rwp', mmm);
            ct_hdr_items    = ct_hdr_items + 1;
        end
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

if ~isnan(opts.MetaDataFieldValues)
    csv_table   = [repmat(opts.MetaDataFieldValues, CakePrms.bins(1), 1) csv_table];
end

% size(csv_table_hdr)
% size(csv_table)

%%% CONVERT TO MATLAB TABLE FORMAT
CSV_TABLE   = array2table(csv_table, 'VariableNames', csv_table_hdr);

if Analysis_Options.save_fits
    disp('###########################')
    fprintf('Saving peak fits in %s\n', pfname_pkfit)
    save(pfname_pkfit, 'fname_pattern', 'Material', 'pkfit', 'pkpars', 'CakePrms', 'Instr', 'polimg');
    
    writetable(CSV_TABLE, pfname_pkfit_tbl)
else
    disp('###########################')
    fprintf('Not saving peak fits for %s\n', fname_pattern)
end
status  = 1;
