function [path_output] = BatchCorrection_RT_function(exp_id, scan_type, ge_num, ...
    root_bkg, bkg_num, ...
    froot, fnum_ini, fnum_fin, fname_metadata, OutAllFrames)
% % clear all
% % close all
% % clc
% 
% % /local/MATLAB/R2019b/bin/matlab -batch "disp(['current folder : ' pwd])"
% 
% % addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% % addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

% OutAllFrames    = true;

%%% START OF THE SCRIPT
path_bkg    = fullfile('/home/beams/S1IDUSER/mnt/s1c', exp_id);
path_image  = fullfile('/home/beams/S1IDUSER/mnt/s1c', exp_id);

% fname_metadata  = 'waxs_hydra_exposures.par';
% pfname_metadata = fullfile(pname_metadata, fname_metadata);
% metadata        = ReadSpecParFile(pfname_metadata, 'Version', 'mpe_standard');

if (scan_type == 5) || (scan_type == 7) || (scan_type == 8)
    if str2double(ge_num) == 1234
        ge_num  = 1:4;
    elseif str2double(ge_num) == 3
        ge_num  = 3;
    elseif str2double(ge_num) == 5
        ge_num  = 5;
    end
    
    path_output = fullfile(sprintf('/home/beams/S1IDUSER/mnt/orthros/%s_bc/', exp_id), froot)
    
    if exist(path_output, 'dir') == 0
        [SUCCESS,~,~]   = mkdir(path_output);
        disp(sprintf('creating %s', path_output))
    end
    
    for iii = 1:1:length(ge_num)
        genum_iii   = ge_num(iii);
        path_bkgi   = fullfile(path_bkg, sprintf('ge%d', ge_num(iii)));
        path_imagei     = fullfile(path_image, sprintf('ge%d', ge_num(iii)));
        path_outputi    = fullfile(path_output, sprintf('ge%d', ge_num(iii)));
        
        if exist(path_outputi, 'dir') == 0
            [SUCCESS,~,~]   = mkdir(path_outputi);
            disp(sprintf('creating %s', path_outputi))
        end
        
        ct  = 0;
        fname_image_ini     = sprintf('%s_%06d.edf.ge%d', froot, fnum_ini(iii), genum_iii);
        fname_image_fin     = sprintf('%s_%06d.edf.ge%d', froot, fnum_fin(iii), genum_iii);
        pfname_image_ini    = fullfile(path_imagei, fname_image_ini);
        pfname_image_fin    = fullfile(path_imagei, fname_image_fin);
        
        %  isfile(pfname_image_ini)
        % isfile(pfname_image_fin)
        % (ct < 600) && ((~isfile(pfname_image_ini)) || (~isfile(pfname_image_fin)))
        while (ct < 600) && ((~isfile(pfname_image_ini)) || (~isfile(pfname_image_fin)))
            ct  = ct + 1;
            disp(sprintf('waiting for file to show up ... %d s', ct))
            disp(sprintf('%s', pfname_image_ini))
            disp(sprintf('%s', pfname_image_fin))
            pause(60)
        end
        
        pname_metadata  = sprintf('/home/beams/S1IDUSER/new_data/%s', exp_id);
        % fname_metadata  = 'saxs_waxs_fmt_fastpar_REALER.par';
        pfname_metadata = fullfile(pname_metadata, fname_metadata);
        metadata        = ReadSpecParFile(pfname_metadata, 'Version', 'saxs_waxs_fmt_fastpar_v6');
        
        idx_froot       = strcmp([metadata.imgprefix], froot);
        idx_fnum_ini    = getfield(metadata, sprintf('det%d_fnum', genum_iii)) == fnum_ini(iii);
        idx_fnum_fin    = getfield(metadata, sprintf('det%d_fnum', genum_iii)) == fnum_fin(iii);
        idx_froot_fnum_ini  = idx_froot & idx_fnum_ini;
        idx_froot_fnum_fin  = idx_froot & idx_fnum_fin;
        idx_isge            = contains(metadata.det_type, 'GE');
        
        idx_froot_fnum_ini  = find(idx_froot_fnum_ini & idx_isge);
        idx_froot_fnum_fin  = find(idx_froot_fnum_fin & idx_isge);
        
        %%% EXTRACT INFORMATION FROM METADATA
        froot_ge_table  = getfield(metadata, sprintf('det%d_fname', genum_iii), ...
            {idx_froot_fnum_ini:idx_froot_fnum_fin});
        det_type_table  = getfield(metadata, 'det_type', ...
            {idx_froot_fnum_ini:idx_froot_fnum_fin});
        fnum_ge_table	= getfield(metadata, sprintf('det%d_fnum', genum_iii), ...
            {idx_froot_fnum_ini:idx_froot_fnum_fin});
        exptime_table	= getfield(metadata, sprintf('det%d_time_per_frame', genum_iii), ...
            {idx_froot_fnum_ini:idx_froot_fnum_fin});
        
        saxs_nframes_table	= metadata.det5_frames_per_file(idx_froot_fnum_ini:idx_froot_fnum_fin);
        saxs_fnum_fin_table = metadata.det5_fnum(idx_froot_fnum_ini:idx_froot_fnum_fin);
        saxs_fnum_ini_table = saxs_fnum_fin_table - saxs_nframes_table + 1;
        saxs_exptime_table	= metadata.det5_time_per_frame(idx_froot_fnum_ini:idx_froot_fnum_fin);
        
        ic_ct_table = [metadata.scaler1_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler2_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler3_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler4_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler5_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler6_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler7_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler8_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler9_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler10_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler11_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler12_val(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            ];
        
        ic_units_table  = [metadata.scaler1_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler2_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler3_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler4_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler5_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler6_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler7_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler8_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler9_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler10_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler11_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.scaler12_units(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            ];
        
        epoch_time_table    = metadata.epoch_time(idx_froot_fnum_ini:idx_froot_fnum_fin);
        
        temp_data_table     = [ ...
            metadata.ev1(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.ev2(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            metadata.ev3(idx_froot_fnum_ini:idx_froot_fnum_fin), ...
            ];
        
        metadata_table  = table(froot_ge_table, det_type_table, fnum_ge_table, exptime_table, ...
            saxs_nframes_table, saxs_fnum_ini_table, saxs_fnum_fin_table, saxs_exptime_table, ...
            ic_units_table, ic_ct_table, epoch_time_table, temp_data_table);
        
        fname_metadata_table    = sprintf('waxs_%s_%d_%d_ge%d_metadata.csv', ...
            froot, fnum_ge_table(1), fnum_ge_table(end), genum_iii);
        pfname_metadata_table   = fullfile(path_output, fname_metadata_table);
        
        disp(sprintf('summary metadata : %s', fname_metadata_table));
        writetable(metadata_table, pfname_metadata_table, 'WriteRowNames', true);
        
        BatchCorrection(path_bkgi, bkg_num(iii), root_bkg, ...
            path_imagei, froot, ...
            genum_iii, path_outputi, ...
            'CorrectAllImages', false, ...
            'OutAllFrames', OutAllFrames, ...
            'DisplayFrames', false, ...
            'lo', fnum_ini(iii), ...
            'hi', fnum_fin(iii), ...
            'NumDigits', 6)
    end
    
    if (scan_type == 8)
        saxs_pname_output   = fullfile(path_output, 'pixirad');
        if exist(saxs_pname_output, 'dir') == 0
            [SUCCESS,~,~]   = mkdir(saxs_pname_output);
        end
        
        for iii = 1:1:length(saxs_nframes_table)
            clear saxs_imdata_out
            saxs_nframes_iii    = saxs_nframes_table(iii);
            saxs_nframes_actual = 0;
            
            % CHECK IF THE FILES EXIST AND FIND ONES THAT EXIST
            for jjj = 1:1:saxs_nframes_iii
                saxs_fnum_jjj    = saxs_fnum_ini_table(iii) + jjj - 1;
                
                saxs_fname_jjj   = sprintf('%s_%06d.tif', froot, saxs_fnum_jjj);
                saxs_pfname_jjj  = fullfile(path_image, 'pixirad', saxs_fname_jjj);
                if isequal(exist(saxs_pfname_jjj, 'file'), 2)
                    saxs_nframes_actual = saxs_nframes_actual + 1;
                    saxs_imdata_jjj     = ReadPixirad(saxs_pfname_jjj, 'version', 'pixi2');
                    
                    if ~exist('imdata_out', 'var')
                        saxs_imdata_out = saxs_imdata_jjj;
                    else
                        saxs_imdata_out = saxs_imdata_out + saxs_imdata_jjj;
                    end
                end
            end
            
            disp(sprintf('*** %s : %d of %d frames', froot, saxs_nframes_actual, saxs_nframes_iii))
            if saxs_nframes_actual > 1
                
%                 saxs_nframes_table	= metadata.det5_frames_per_file(idx_froot_fnum_ini:idx_froot_fnum_fin);
%                 saxs_fnum_fin_table = metadata.det5_fnum(idx_froot_fnum_ini:idx_froot_fnum_fin);
%                 saxs_fnum_ini_table = saxs_fnum_fin_table - saxs_nframes_table + 1;
%                 saxs_exptime_table	= metadata.det5_time_per_frame(idx_froot_fnum_ini:idx_froot_fnum_fin);
%                 
                saxs_imdata_out_ave     = saxs_imdata_out./saxs_nframes_actual;
                saxs_imdata_out_per_sec = saxs_imdata_out_ave/saxs_exptime_table(iii);
                
                %%% TIFF FILE OUT
                saxs_pname_output_tif   = fullfile(saxs_pname_output, 'ave');
                if exist(saxs_pname_output_tif, 'dir') == 0
                    [SUCCESS,~,~]   = mkdir(saxs_pname_output_tif);
                end
                saxs_fname_output_tif   = sprintf('%s_%06d_%06d_ave.tif', ...
                    froot, saxs_fnum_ini_table(iii), saxs_fnum_fin_table(iii));
                saxs_pfname_output_tif  = fullfile(saxs_pname_output_tif, saxs_fname_output_tif);
                %%% IGOR / NIKA WANTS NO FLIP
                WritePixiradAveTIFF(saxs_pfname_output_tif, saxs_imdata_out_ave, ...
                    'apply_flip', false)
                
                %%% TIFF FILE OUT
                saxs_pname_output_tif    = fullfile(saxs_pname_output ,'ave_1s');
                if exist(saxs_pname_output_tif, 'dir') == 0
                    [SUCCESS,~,~]   = mkdir(saxs_pname_output_tif);
                end
                saxs_fname_output_tif    = sprintf('%s_%06d_%06d_ave_1s_norm.tif', ...
                    froot, saxs_fnum_ini_table(iii), saxs_fnum_fin_table(iii));
                saxs_pfname_output_tif   = fullfile(saxs_pname_output_tif, saxs_fname_output_tif);
                %%% IGOR / NIKA WANTS NO FLIP
                WritePixiradAveTIFF(saxs_pfname_output_tif, saxs_imdata_out_per_sec, ...
                    'apply_flip', false)
                
                %%% MAT FILE OUT
                saxs_pname_output_mat   = fullfile(saxs_pname_output ,'mat');
                if exist(saxs_pname_output_mat, 'dir') == 0
                    [SUCCESS,~,~]   = mkdir(saxs_pname_output_mat);
                end
                saxs_fname_output_mat   = sprintf('%s_%06d_%06d_ave_1s_norm.mat', ...
                    froot, saxs_fnum_ini_table(iii), saxs_fnum_fin_table(iii));
                saxs_pfname_output      = fullfile(saxs_pname_output_mat, saxs_fname_output_mat);
                save(saxs_pfname_output, 'saxs_imdata_out_ave', 'saxs_imdata_out_per_sec');
                
                disp(sprintf('output files take the form of %s', saxs_fname_output_mat));
                
                saxs_fname_output_tif_table{iii, 1} = saxs_fname_output_tif;
                
                clear imdata_out
            else
                saxs_fname_output_tif_table{iii, 1}  = 'NIL';
            end
        end
        
        saxs_metadata_table = table(froot_ge_table, saxs_fname_output_tif_table, fnum_ge_table, exptime_table, ...
            saxs_nframes_table, saxs_fnum_ini_table, saxs_fnum_fin_table, saxs_exptime_table, ...
            ic_units_table, ic_ct_table, epoch_time_table, temp_data_table);
        saxs_fname_metadata_table   = sprintf('saxs_%s_%d_%d_ge_%d_%d_pixirad_metadata.csv', ...
            froot, fnum_ge_table(1), fnum_ge_table(end), ...
            saxs_fnum_ini_table(1), saxs_fnum_fin_table(end));
        saxs_pfname_metadata_table  = fullfile(path_output, saxs_fname_metadata_table);
        
        disp(sprintf('summary metadata : %s', saxs_fname_metadata_table));
        writetable(saxs_metadata_table, saxs_pfname_metadata_table, 'WriteRowNames', true);
    end
end