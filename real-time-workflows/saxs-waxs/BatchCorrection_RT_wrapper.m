clear all
close all
clc

addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%% BC OUTPUT ALL FRAMES
OutAllFrames    = false;

%%% BC OUTPUT ALL FRAMES
DoLineout       = true;


%%% PVs check
pv_exp_id       = '1ide1:userStringCalc10.KK';
pv_ge_num       = '1ide1:userStringCalc10.LL';

pv_scan_type    = '1ide1:userStringCalc10.K';
pv_scan_state   = '1ide1:userStringCalc10.L';
pv_root_bkg     = '1ide1:userStringCalc10.AA';
pv_bkg_num      = '1ide1:userStringCalc10.BB';

pv_froot        = '1ide1:userStringCalc10.CC';
pv_fnum_ini     = '1ide1:userStringCalc10.DD';
pv_fnum_fin     = '1ide1:userStringCalc10.EE';

%%%
cmdstr      = sprintf('caget -t %s', pv_exp_id);
[~, exp_id] = unix(cmdstr);
exp_id      = strtrim(exp_id);

cmdstr      = sprintf('caget -t %s', pv_ge_num);
[~, ge_num] = unix(cmdstr);
ge_num      = strtrim(ge_num);

bc_job_name_list = {};
ct_job  = 1;
while true 
    disp('*** checking scan state');
    cmdstr          = sprintf('caget -t %s', pv_scan_state);
    [~, scan_state] = unix(cmdstr);
    scan_state      = str2double(scan_state);
    
    cmdstr          = sprintf('caget -t %s', pv_scan_type);
    [~, scan_type] = unix(cmdstr);
    scan_type       = str2double(scan_type);
    
    if (scan_state == 606060)
        disp('****** launching background correction if a new job')
        
        cmdstr          = sprintf('caget -t %s', pv_bkg_num);
        [~, bkg_num]    = unix(cmdstr);
        bkg_num         = str2num(bkg_num);
        
        cmdstr          = sprintf('caget -t %s', pv_root_bkg);
        [~, root_bkg]   = unix(cmdstr);
        root_bkg        = strtrim(root_bkg);
        
        cmdstr      = sprintf('caget -t %s', pv_froot);
        [~, froot]  = unix(cmdstr);
        froot       = strtrim(froot);
        
        cmdstr          = sprintf('caget -t %s', pv_fnum_ini);
        [~, fnum_ini]   = unix(cmdstr);
        fnum_ini        = str2num(fnum_ini);
        
        cmdstr          = sprintf('caget -t %s', pv_fnum_fin);
        [~, fnum_fin]   = unix(cmdstr);
        fnum_fin        = str2num(fnum_fin);
        
        if fnum_fin >= fnum_ini
            bc_job_name = sprintf('%s.%s.fini=%d.ffin=%d.ge%s', ...
                exp_id, froot, fnum_ini(1), fnum_fin(1), ge_num);

            if sum(strcmp(bc_job_name_list, bc_job_name)) == 0
                disp(sprintf('%s is a new job', bc_job_name))
                
                bc_job_name_list{ct_job,1}  = bc_job_name;
                
                job(ct_job) = batch(@BatchCorrection_RT_function, 1, ...
                    {exp_id, scan_type, ge_num, root_bkg, bkg_num, froot, fnum_ini, fnum_fin, ...
                    'saxs_waxs_fmt_fastpar_REALER.par', OutAllFrames});

                ct_job  = ct_job + 1;
                
                if DoLineout
                    disp(sprintf('generate lineout for %s', bc_job_name))
                end
            else
                disp(sprintf('%s was run already', bc_job_name))
            end
        else
            disp('scan captured in the middle')
        end
    elseif (scan_state == 6060)
        disp('****** layer scan on going')
        %%% CLEAN UP SOME JOBS
%         if exist('job', 'var')
%             idx_finished    = false(length(job),1);
%             idx_deleted     = false(length(job),1);
%             
%             for iii = 1:1:length(job)
%                 idx_finished(iii,1) = strcmp(job(iii).State, 'finished');
%                 idx_deleted(iii,1)     = strcmp({job(:).State}, 'deleted');
%                 idx_to_delete   = idx_finished & ~idx_deleted;
%             end
%             
%             if sum(idx_to_delete) > 0
%                 disp('********* deleting finished jobs to free up resources')
%                 job(idx_to_delete).delete;
%             end
%         end
    else
        disp('****** scan doing something else')
    end
    pause(0.2)
end
return

pv_bkg_num  = '1ide1:userStringCalc10.A';
cmdstr      = sprintf('caget -t %s', pv_bkg_num);

[~, bkg_num]    = unix(cmdstr);
bkg_num         = str2double(bkg_num);

pv_root_bkg = '1ide1:userStringCalc10.AA';
cmdstr      = sprintf('caget -t %s', pv_root_bkg);
[~, root_bkg]   = unix(cmdstr);
root_bkg        = strtrim(root_bkg);

return
while true
    
end
path_bkg    = '/home/beams/S1IDUSER/mnt/s1c/connolly_oct22';
path_image  = '/home/beams/S1IDUSER/mnt/s1c/connolly_oct22';

pname_metadata  = '/home/beams/S1IDUSER/new_data/connolly_oct22';
% pname_metadata  = 'D:\s\connolly_oct22\metadata\connolly_oct22';

fname_metadata  = 'saxs_waxs_fmt_fastpar_new.par';
pfname_metadata = fullfile(pname_metadata, fname_metadata);
metadata        = ReadSpecParFile(pfname_metadata, 'Version', 'saxs_waxs_fmt_fastpar_v6');

% fname_metadata  = 'waxs_hydra_exposures.par';
% pfname_metadata = fullfile(pname_metadata, fname_metadata);
% metadata        = ReadSpecParFile(pfname_metadata, 'Version', 'mpe_standard');


root_image
return



genum       = 1:4;   %%% GET FROM A PV IN USER STRCALC10 
bkg_num     = 9426;  %%%f8547
root_bkg    = 'dark_before';
root_image      = 'EB_5_300cycles'; %%%%
OutAllFrames    = true;

% /local/MATLAB/R2019b/bin/matlab -batch "disp(['current folder : ' pwd])"











%%%
% path_output = fullfile('/home/beams/S1IDUSER/mnt/orthros/connolly_oct22_bc/', root_image);
path_output = fullfile('/home/beams/S1IDUSER/mnt/s1c/connolly_oct22_bc/', root_image);
%%%%

if exist(path_output, 'dir') == 0
    [SUCCESS,~,~]   = mkdir(path_output)
end
    
%%%% FILTER OUT METADATA
idx = find(strncmp([metadata.imgprefix], root_image, length(root_image)));

istomo      = false;
ispixirad   = false;

for jjj = 1:1:length(idx)
    idx_jjj     = idx(jjj);
    
    %%% USING SAXS_WAXS UNIFIED PAR FILE
    det_type_jjj    = metadata.det_type{idx_jjj};
    if strcmp(det_type_jjj, 'PointGrey')
        istomo  = true;
    elseif strcmp(det_type_jjj, 'GE_STILL') || strcmp(det_type_jjj, 'GE_AD') || strcmp(det_type_jjj, 'GE_SW_FS')
        
        %%% USING SAXS_WAXS UNIFIED PAR FILE
        fnum_ge1_jjj    = metadata.det1_fnum(idx_jjj);
        exptime1_jjj    = metadata.det1_time_per_frame(idx_jjj);
        fnum_ge2_jjj    = metadata.det2_fnum(idx_jjj);
        exptime2_jjj    = metadata.det2_time_per_frame(idx_jjj);
        fnum_ge3_jjj    = metadata.det3_fnum(idx_jjj);
        exptime3_jjj    = metadata.det3_time_per_frame(idx_jjj);
        fnum_ge4_jjj    = metadata.det4_fnum(idx_jjj);
        exptime4_jjj    = metadata.det4_time_per_frame(idx_jjj);
        fnum_ge5_jjj    = metadata.det7_fnum(idx_jjj);
        exptime5_jjj    = metadata.det7_time_per_frame(idx_jjj);
        
        saxs_nframes_jjj    = metadata.det5_frames_per_file(idx_jjj);
        fnum_saxs_fin       = metadata.det5_fnum(idx_jjj);
        fnum_saxs_ini       = fnum_saxs_fin - saxs_nframes_jjj + 1;
        exptime_saxs_jjj    = metadata.det5_time_per_frame(idx_jjj);
        
        %%% WAXS_SAXS PAR
        % fnum_fin    = metadata.det5_frames_per_file(idx_jjj);
        % fnum_ini    = metadata.det5_fnum(idx_jjj);
        % nframes_jjj = fnum_fin - fnum_ini + 1;
        %%%
        
        %%% ASSEMBLE METADATA FOR A DATA POINT
        %         fnum_jjj    = [fnum_ge1_jjj, ...
        %             fnum_ge2_jjj, ...
        %             fnum_ge3_jjj, ...
        %             fnum_ge4_jjj, ...
        %             fnum_ge5_jjj, ...
        %             fnum_saxs_ini, ...
        %             ];
        
        fnum_jjj    = [fnum_ge1_jjj, ...
            fnum_ge5_jjj, ...
            fnum_saxs_ini, ...
            ];
        
        % exptime_jjj = [exptime1_jjj exptime2_jjj exptime3_jjj exptime4_jjj exptime5_jjj exptime_saxs_jjj];
        exptime_jjj = [exptime1_jjj exptime5_jjj exptime_saxs_jjj];
        
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
        
        %%% COMPILE THE LIST
        ic_ct_list(jjj, :)      = ic_ct_jjj;
        ic_units_list(jjj, :)   = ic_units_jjj;
        exptime_list(jjj, :)    = exptime_jjj;
        epoch_time(jjj, 1)      = epoch_time_jjj;
        temperature_data(jjj,:) = temperature_data_jjj;
        fnum_list(jjj,:)        = fnum_jjj;
        
        det_type{jjj,1}             = det_type_jjj;
        froot_output_list{jjj, 1}   = root_image;
    end
end

if ~istomo
    tableout    = table(froot_output_list, det_type, fnum_list, exptime_list, ...
        ic_units_list, ic_ct_list, epoch_time, temperature_data);
    
    fname_table     = sprintf('%s_saxs_waxs_metadata.csv', root_image);
    pfname_table    = fullfile(path_output, fname_table)
    
    disp(sprintf('summary metadata : %s', fname_table));
    writetable(tableout, pfname_table, 'WriteRowNames', true);
end

for iii = 1:1:length(genum)
    genum_iii   = genum(iii);
    path_bkgi   = fullfile(path_bkg, sprintf('ge%d', genum(iii)));
    path_imagei     = fullfile(path_image, sprintf('ge%d', genum(iii)));
    path_outputi    = fullfile(path_output, sprintf('ge%d', genum(iii)));
    
    if exist(path_outputi, 'dir') == 0
        [SUCCESS,~,~]   = mkdir(path_outputi)
    end
    
    BatchCorrection(path_bkgi, bkg_num, root_bkg, ...
        path_imagei, root_image, ...
        genum_iii, path_outputi, ...
        'CorrectAllImages', true, ...
        'OutAllFrames', OutAllFrames, ...
        'DisplayFrames', false, ...
        'NumDigits', 6)
end