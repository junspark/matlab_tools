clear all
close all
clc

addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

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
XRDIMAGE.ExpID              = 'startup_feb22';
XRDIMAGE.MetaDataFile       = '/home/beams/S1IDUSER/new_data/startup_feb22/waxs_hydra_exposures.par';
XRDIMAGE.Image.pname_chi    = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/startup_feb22_analysis/waxs/stock';

XRDIMAGE.Image.pname        = '/net/wolf/data/tomo1/startup_feb22_bc/stock_dogfish_scan1';
XRDIMAGE.Image.fbase        = 'stock_dogfish_scan1';
XRDIMAGE.Image.fnumber      = 157:247;

% XRDIMAGE.Image.pname        = '/net/wolf/data/tomo1/startup_feb22_bc/stock_dogfish_scan2';
% XRDIMAGE.Image.fbase        = 'stock_dogfish_scan2';
% XRDIMAGE.Image.fnumber      = 250:340;

% XRDIMAGE.Image.pname        = '/net/wolf/data/tomo1/startup_feb22_bc/stock_dogfish_scan3';
% XRDIMAGE.Image.fbase        = 'stock_dogfish_scan3';
% XRDIMAGE.Image.fnumber      = 343:358;

XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.corrected    = 1;        % 0 - UNCORRECTED, 1 = SUM, 2 = COR32
XRDIMAGE.Image.dettype      = 5;        % 3 = GE3, 5 = GE5, 1234 = HYDRA1234

%%% ANALYSIS STARTS HERE
metadata    = ReadSpecParFile(XRDIMAGE.MetaDataFile, 'Version', 'mpe_standard');

switch XRDIMAGE.Image.dettype
    case 1234
        sprintf('%d = HYDRA1234', XRDIMAGE.Image.dettype);
        
        for iii = 1:1:4
            fext    = sprintf('ge%d', iii);
            fimctrl = sprintf('det_%s_72keV_2030mm_caking.imctrl', fext);
            
            imctrl{iii} = ReadGSAS2Imctrl(fimctrl);
            
            azm_range   = imctrl{iii}.azm_range;
            azm_offset  = imctrl{iii}.azm_offset;
            azm_steps   = imctrl{iii}.azm_steps;
            tth_range   = imctrl{iii}.tth_range;
            tth_steps   = imctrl{iii}.tth_steps;
            
            dazm        = (azm_range(2) - azm_range(1))/azm_steps;
            switch iii
                case 1
                    azm_grid    = linspace(azm_range(1)+dazm/2, azm_range(2)-dazm/2, azm_steps) + 67;
                otherwise
                    azm_grid    = linspace(azm_range(1)+dazm/2, azm_range(2)-dazm/2, azm_steps);
            end
            
            idx_azm                 = azm_grid >= 360;
            azm_grid(idx_azm)       = azm_grid(idx_azm) - 360;
            imctrl{iii}.azm_grid    = azm_grid;
        end
        
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            
            idx_fbase       = ismember(metadata.det1_fname, XRDIMAGE.Image.fbase);
            idx_fnumber     = ismember(metadata.det1_fnum, XRDIMAGE.Image.fnumber(iii));
            idx_metadata    = find(idx_fbase & idx_fnumber);
            switch ~isempty(idx_metadata) && (length(idx_metadata) == 1)
                case true
                    nframes         = metadata.det1_frames_per_file(idx_metadata);
                    
                    %%% PARFOR HERE
                    parfor jjj = 0:nframes
                        ydata_sum   = zeros(imctrl{1}.tth_steps, 1);
                        %%% LOOP OVER GE 1 ... 4
                        for kkk = 1:1:4
                            fext    = sprintf('ge%d', kkk);
                            
                            for mmm = 1:1:imctrl{kkk}.azm_steps
                                switch jjj
                                    case 0
                                        fchi    = sprintf('%s.%s.sum_Azm=_%d.chi', froot, fext, fix(imctrl{kkk}.azm_grid(mmm)));
                                    otherwise
                                        fchi    = sprintf('%s.%s_frame_%d.cor32_Azm=_%d.chi', froot, fext, jjj, fix(imctrl{kkk}.azm_grid(mmm)));
                                end
                                pfchi   = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fext, fchi);
                                
                                switch isfile(pfchi)
                                    case true
                                        disp(sprintf('%s DOES EXISTS\n', fchi));
                                        pfchi_data_kkk_mmm  = table2array(readtable(pfchi, 'FileType', 'text', ...
                                            'ReadVariableNames', false, ...
                                            'HeaderLines', 4));
                                        ydata_sum   = ydata_sum + pfchi_data_kkk_mmm(:,2);
                                    otherwise
                                        disp(sprintf('missing %s, re-run auto-integration', fchi))
                                end
                                  
                            end
                        end
                        xdata       = pfchi_data_kkk_mmm(:,1);
                        tableout    = [xdata ydata_sum];
                        switch jjj
                            case 0
                                fchi_out    = sprintf('%s.sum.chi', froot);
                            otherwise
                                fchi_out    = sprintf('%s.frame_%d.cor32.chi', froot, jjj);
                        end
                        pfchi_out   = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fchi_out);
                        
                        fid = fopen(pfchi_out, 'w');
                        fprintf(fid, '%s Azm= 0.0\n', fchi_out);
                        fprintf(fid, '2-Theta Angle (Degrees)\n');
                        fprintf(fid, 'Intensity\n');
                        fprintf(fid, '       %d\n', imctrl{1}.tth_steps);
                        fclose(fid);
                        dlmwrite(pfchi_out, tableout, '-append', ...
                            'precision', '%.7f', ...
                            'newline', 'unix', ...
                            'delimiter', ' ', ...
                            'coffset', 1)
                        % plot(xdata, ydata, pfchi_data_kkk_mmm(:,1), pfchi_data_kkk_mmm(:,2));
                    end
                otherwise
                    warning('%s metadata needs investigation')
            end
        end
    case 5
        for iii = 5
            fext    = sprintf('ge%d', iii);
            fimctrl = 'startup_feb22_71keV_800mm_stock.imctrl';
            
            imctrl{iii} = ReadGSAS2Imctrl(fimctrl);
            
            azm_range   = imctrl{iii}.azm_range;
            azm_offset  = imctrl{iii}.azm_offset;
            azm_steps   = imctrl{iii}.azm_steps;
            tth_range   = imctrl{iii}.tth_range;
            tth_steps   = imctrl{iii}.tth_steps;
            
            dazm        = (azm_range(2) - azm_range(1))/azm_steps;
            switch iii
                case 1
                    azm_grid    = linspace(azm_range(1)+dazm/2, azm_range(2)-dazm/2, azm_steps) + 67;
                otherwise
                    azm_grid    = linspace(azm_range(1)+dazm/2, azm_range(2)-dazm/2, azm_steps);
            end
            
            idx_azm                 = azm_grid >= 360;
            azm_grid(idx_azm)       = azm_grid(idx_azm) - 360;
            imctrl{iii}.azm_grid    = azm_grid;
        end
        
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            
            idx_fbase       = ismember(metadata.det6_fname, XRDIMAGE.Image.fbase);
            idx_fnumber     = ismember(metadata.det6_fnum, XRDIMAGE.Image.fnumber(iii));
            idx_metadata    = find(idx_fbase & idx_fnumber);
            switch ~isempty(idx_metadata) && (length(idx_metadata) == 1)
                case true
                    nframes         = metadata.det1_frames_per_file(idx_metadata);
                    
                    %%% PARFOR HERE
                    % for jjj = 0:nframes
                    parfor jjj = 0%:nframes
                        ydata_sum   = zeros(imctrl{5}.tth_steps, 1);
                        ydata_azm   = [];
                        %%% LOOP OVER GE 1 ... 4
                        for kkk = 5
                            fext    = sprintf('ge%d', kkk);
                            for mmm = 1:1:imctrl{kkk}.azm_steps
                                switch jjj
                                    case 0
                                        fchi    = sprintf('%s.edf.%s.sum_Azm=_%d.chi', froot, fext, fix(imctrl{kkk}.azm_grid(mmm)));
                                    otherwise
                                        fchi    = sprintf('%s.%s_frame_%d.cor32_Azm=_%d.chi', froot, fext, jjj, fix(imctrl{kkk}.azm_grid(mmm)));
                                end
                                pfchi   = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fext, fchi);
                                
                                switch isfile(pfchi)
                                    case true
                                        disp(sprintf('%s DOES EXISTS\n', fchi));
                                        pfchi_data_kkk_mmm  = table2array(readtable(pfchi, 'FileType', 'text', ...
                                            'ReadVariableNames', false, ...
                                            'HeaderLines', 4));
                                        ydata_sum   = ydata_sum + pfchi_data_kkk_mmm(:,2);
                                        ydata_azm   = [ydata_azm pfchi_data_kkk_mmm(:,2)];
                                    otherwise
                                        disp(sprintf('missing %s, re-run auto-integration', fchi))
                                end
                            end
                        end
                        
                        xdata               = pfchi_data_kkk_mmm(:,1);
                        tableout_sum_azm    = [xdata ydata_sum];
                        tableout_indiv_azm  = [xdata ydata_azm];
                        switch jjj
                            case 0
                                fchi_out_sum_azm    = sprintf('%s.sum_azm.chi', froot);
                                fchi_out_ind_azm    = sprintf('%s.ind_azm.chi', froot);
                            otherwise
                                fchi_out_sum_azm    = sprintf('%s.frame_%d.cor32.sum_azm.chi', froot, jjj);
                                fchi_out_ind_azm    = sprintf('%s.frame_%d.cor32.ind_azm.chi', froot, jjj);
                        end
                        pfchi_out_sum_azm   = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fchi_out_sum_azm);
                        pfchi_out_ind_azm   = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fchi_out_ind_azm);
                        
                        fid = fopen(pfchi_out_sum_azm, 'w');
                        fprintf(fid, '%s Azm= 0.0\n', fchi_out_sum_azm);
                        fprintf(fid, '2-Theta Angle (Degrees)\n');
                        fprintf(fid, 'Intensity\n');
                        fprintf(fid, '       %d\n', imctrl{5}.tth_steps);
                        fclose(fid);
                        dlmwrite(pfchi_out_sum_azm, tableout_sum_azm, '-append', ...
                            'precision', '%.7f', ...
                            'newline', 'unix', ...
                            'delimiter', ' ', ...
                            'coffset', 1)
                        
                        fid = fopen(pfchi_out_ind_azm, 'w');
                        fprintf(fid, '%s Azm= ', fchi_out_sum_azm);
                        fprintf(fid, '%f ', imctrl{kkk}.azm_grid);
                        fprintf(fid, '\n');
                        fprintf(fid, '2-Theta Angle (Degrees)\n');
                        fprintf(fid, 'Intensity\n');
                        fprintf(fid, '       %d\n', imctrl{5}.tth_steps);
                        fclose(fid);
                        dlmwrite(pfchi_out_ind_azm, tableout_indiv_azm, '-append', ...
                            'precision', '%.7f', ...
                            'newline', 'unix', ...
                            'delimiter', ' ', ...
                            'coffset', 1)
                            
                        % plot(xdata, ydata, pfchi_data_kkk_mmm(:,1), pfchi_data_kkk_mmm(:,2));
                    end
                otherwise
                    warning('%s metadata needs investigation')
            end
        end
    otherwise
        warning('%d is not a valid option', XRDIMAGE.Image.dettype);
end
