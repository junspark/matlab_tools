clear all
close all
clc

% addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

% path_bkg    = '/home/beams/S1IDUSER/mnt/s1c/mpe_jun22/pilatus';
% path_image  = '/home/beams/S1IDUSER/mnt/s1c/mpe_jun22/pilatus';
path_bkg    = 'D:\s\mpe_jun22\pilatus';
path_image  = 'D:\s\mpe_jun22\pilatus';

% pname_metadata  = '/home/beams/S1IDUSER/new_data/mpe_jun22';
pname_metadata  = 'D:\s\mpe_jun22\metadata\mpe_jun22';

fname_metadata  = 'saxs_waxs_fmt_fastpar_new.par';
pfname_metadata = fullfile(pname_metadata, fname_metadata);
metadata        = ReadSpecParFile(pfname_metadata, 'Version', 'saxs_waxs_fmt_fastpar_v6');

% fname_metadata  = 'waxs_hydra_exposures.par';
% pfname_metadata = fullfile(pname_metadata, fname_metadata);
% metadata        = ReadSpecParFile(pfname_metadata, 'Version', 'mpe_standard');

% pname_fout      = '/home/beams/S1IDUSER/mnt/orthros/mpe_jun22_bc/';
pname_fout      = 'D:\orthros\mpe_jun22_bc';

froot_image   = unique([metadata.imgprefix]);

for iii = 15:1:length(froot_image)
    froot_image_iii = froot_image{iii};
    disp(sprintf('******** %s ********', froot_image_iii));
    
    %%% USING SAXS_WAXS UNIFIED PAR FILE
    idx = find(strcmp([metadata.imgprefix], froot_image_iii));
    %%% WAXS_SAXS PAR
    % idx = find(strcmp([metadata.det5_fname], froot_image_iii));
    
    %     ic_ct_list      = [];
    %     ic_units_list   = [];
    %     pfname_output_tif_list  = [];
    clear ic_ct_list ic_units_list exptime epoch_time temperature_data fname_output_tif_list
                
    ispilatus   = false;
    
    for jjj = 1:1:length(idx)
        idx_jjj     = idx(jjj);
        
        %%% USING SAXS_WAXS UNIFIED PAR FILE
        det_type    = metadata.det_type{idx_jjj};
        
        if contains(det_type, 'pilatus', 'IgnoreCase', true)
            ispilatus   = true;
            
            disp(sprintf('%s, %s', froot_image_iii, det_type))
            
            path_output_iii = fullfile(pname_fout, froot_image_iii, 'pilatus');
            if exist(path_output_iii, 'dir') == 0
                [SUCCESS,~,~]   = mkdir(path_output_iii);
            end
            
            %%% USING SAXS_WAXS UNIFIED PAR FILE
            nframes_jjj = metadata.det8_fnum(idx_jjj);
            fnum_fin    = metadata.det8_frames_per_file(idx_jjj);
            fnum_ini    = fnum_fin - nframes_jjj + 1;
            exptime_jjj = metadata.det8_time_per_frame(idx_jjj);
            
            %%% WAXS_SAXS PAR
            % fnum_fin    = metadata.det5_frames_per_file(idx_jjj);
            % fnum_ini    = metadata.det5_fnum(idx_jjj);
            % nframes_jjj = fnum_fin - fnum_ini + 1;
            
            
            disp(sprintf('fini = %d, fend = %d, nframes = %d, exp time = %g sec', fnum_ini, fnum_fin, nframes_jjj, exptime_jjj))
            
            ic_ct_jjj   = [metadata.scaler1_val(idx_jjj), ...
                metadata.scaler2_val(idx_jjj), ...
                metadata.scaler3_val(idx_jjj), ...
                metadata.scaler4_val(idx_jjj), ...
                metadata.scaler5_val(idx_jjj), ...
                metadata.scaler6_val(idx_jjj), ...
                metadata.scaler7_val(idx_jjj), ...
                metadata.scaler8_val(idx_jjj), ...
                metadata.scaler9_val(idx_jjj), ...
                metadata.scaler10_val(idx_jjj), ...
                metadata.scaler11_val(idx_jjj), ...
                metadata.scaler12_val(idx_jjj), ...
                ];
            
            ic_units_jjj    = [metadata.scaler1_units(idx_jjj), ...
                metadata.scaler2_units(idx_jjj), ...
                metadata.scaler3_units(idx_jjj), ...
                metadata.scaler4_units(idx_jjj), ...
                metadata.scaler5_units(idx_jjj), ...
                metadata.scaler6_units(idx_jjj), ...
                metadata.scaler7_units(idx_jjj), ...
                metadata.scaler8_units(idx_jjj), ...
                metadata.scaler9_units(idx_jjj), ...
                metadata.scaler10_units(idx_jjj), ...
                metadata.scaler11_units(idx_jjj), ...
                metadata.scaler12_units(idx_jjj), ...
                ];
            epoch_time_jjj          = metadata.epoch_time(idx_jjj);
            temperature_data_jjj    = [ ...
                metadata.ev1(idx_jjj), ...
                metadata.ev2(idx_jjj), ...
                metadata.ev3(idx_jjj), ...
                ];
            
            nframes_actual	= 0;
            for kkk = 1:1:nframes_jjj
                fnum_kkk    = fnum_ini + kkk - 1;
                fname_kkk   = sprintf('%s_%06d.tif', froot_image_iii, fnum_kkk);
                pfname_kkk  = fullfile(path_image, froot_image_iii, fname_kkk);
                if isequal(exist(pfname_kkk, 'file'), 2)
                    nframes_actual  = nframes_actual + 1;
                    imdata_kkk      = ReadPilatus(pfname_kkk);
                    
                    % imagesc(imdata_kkk)
                    % axis equal tight
                    % caxis([0 10])
                    if ~exist('imdata_out', 'var')
                        imdata_out  = imdata_kkk;
                    else
                        imdata_out  = imdata_out + imdata_kkk;
                    end
                end
            end
            
            disp(sprintf('%s : %d frames', froot_image_iii, nframes_actual))
            if nframes_actual > 1
                imdata_out_ave      = imdata_out./nframes_actual;
                imdata_out_per_sec  = imdata_out_ave/exptime_jjj;
                
                %%% TIFF FILE OUT
                pname_output_tif_iii    = fullfile(path_output_iii ,'ave');
                if exist(pname_output_tif_iii, 'dir') == 0
                    [SUCCESS,~,~]   = mkdir(pname_output_tif_iii);
                end
                fname_output_tif_iii    = sprintf('%s_%06d_%06d_ave.tif', ...
                    froot_image_iii, fnum_ini, fnum_fin);
                pfname_output_tif_iii   = fullfile(pname_output_tif_iii, fname_output_tif_iii);
                
                WritePilatusAveTIFF(pfname_output_tif_iii, imdata_out_ave)
                
                %%% TIFF FILE OUT
                pname_output_tif_iii    = fullfile(path_output_iii ,'ave_1s');
                if exist(pname_output_tif_iii, 'dir') == 0
                    [SUCCESS,~,~]   = mkdir(pname_output_tif_iii);
                end
                fname_output_tif_iii    = sprintf('%s_%06d_%06d_ave_1s_norm.tif', ...
                    froot_image_iii, fnum_ini, fnum_fin);
                pfname_output_tif_iii   = fullfile(pname_output_tif_iii, fname_output_tif_iii);
                
                WritePilatusAveTIFF(pfname_output_tif_iii, imdata_out_per_sec)
                
                %%% MAT FILE OUT
                pname_output_mat_iii    = fullfile(path_output_iii ,'mat');
                if exist(pname_output_mat_iii, 'dir') == 0
                    [SUCCESS,~,~]   = mkdir(pname_output_mat_iii);
                end
                fname_output_mat_iii    = sprintf('%s_%06d_%06d_ave_1s_norm.mat', ...
                    froot_image_iii, fnum_ini, fnum_fin);
                pfname_output_iii   = fullfile(pname_output_mat_iii, fname_output_mat_iii);
                save(pfname_output_iii, 'imdata_out_ave', 'imdata_out_per_sec');
                
                disp(sprintf('output files take the form of %s', fname_output_mat_iii));
                
                ic_ct_list(jjj, :)      = ic_ct_jjj;
                ic_units_list(jjj, :)   = ic_units_jjj;
                exptime(jjj, 1)         = exptime_jjj;
                epoch_time(jjj, 1)      = epoch_time_jjj;
                temperature_data(jjj,:) = temperature_data_jjj;
                
                fname_output_tif_list{jjj, 1}  = fname_output_tif_iii;
                
                clear imdata_out
            else
                ic_ct_list(jjj, :)      = ic_ct_jjj*nan;
                ic_units_list(jjj, :)   = ic_units_jjj*nan;
                exptime(jjj, 1)         = nan;
                epoch_time(jjj, 1)      = nan;
                temperature_data(jjj,:) = temperature_data_jjj*nan;
                fname_output_tif_list{jjj, 1}  = 'NIL';
            end
        else
            disp(sprintf('%s, %s', froot_image_iii, det_type))
            disp('skipping ... ')
        end
    end
    
    if ispilatus
        tableout    = table(fname_output_tif_list, epoch_time, exptime, ic_units_list, ic_ct_list, temperature_data);
        
        fname_table     = sprintf('%s_rapid_access_metadata.csv', froot_image_iii);
        pfname_table    = fullfile(path_output_iii, fname_table);
        
        disp(sprintf('summary metadata : %s', fname_table));
        writetable(tableout, pfname_table, 'WriteRowNames', true);
    end
end
