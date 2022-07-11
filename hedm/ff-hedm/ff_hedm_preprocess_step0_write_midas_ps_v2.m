%%%% THIS IS BASED ON KASEMER_JUN22 DATA

clear all
close all
clc

%%%%%
addpath(genpath('/home/beams0/PARKJS/matlab/matlab_tools'));

pname_data      = '/home/beams/S1IDUSER/mnt/s1c/kasemer_jun22/ge5';
data_fext       = 'edf.ge5';
data_ome_step_size  = 0.25;

pname_metadata  = '/home/beams/S1IDUSER/new_data/kasemer_jun22';
fname_metadata  = 'saxs_waxs_fmt_fastpar_new.par';
pfname_metadata = fullfile(pname_metadata, fname_metadata);
metadata        = ReadSpecParFile(pfname_metadata, 'Version', 'saxs_waxs_fmt_fastpar_v6');

pname_ps_template   = '/home/beams/S1IDUSER/mnt/orthros/kasemer_jun22_midas/ff/shankar/';
fname_ps_template   = 'ps_Ti_template.txt';
pname_ps_target     = pname_ps_template;

%%%%%
froot_imageprefix   = unique([metadata.imgprefix]);
froot_sample_name   = 'sam1';
% froot_sample_name   = 'sam3';

fname_midas_script  = sprintf('midas_run.%s.out', froot_sample_name);

ct  = 0;
for iii = 11:1:length(froot_imageprefix)
    froot_image_iii = froot_imageprefix{iii};
    disp(sprintf('******** %s ********', froot_image_iii));
    
    idx = find(strcmp([metadata.imgprefix], froot_image_iii));
    
    %%% CHECK IF THE SCAN IS FF SCAN USING GE
    isge    = true;
    for jjj = 1:1:length(idx)
        if ~strcmp(metadata.det_type{idx(jjj)}, 'GE_AD')
            isge    = false;
        end
    end
    
    %%% CHECK IF THE GE FILES EXIST AND CONSISTENT WITH METADATA
    is_good_ff_scan = true;
    if isge
        for jjj = 1:1:length(idx)
            idx_jjj     = idx(jjj);
            fnum_jjj    = metadata.det7_fnum(idx_jjj);
            
            fname_gefile    = sprintf('%s_%06d.%s', froot_image_iii, fnum_jjj, data_fext);
            pfname_gefile   = fullfile(pname_data, fname_gefile);
            
            nframes_per_file    = metadata.det7_frames_per_file(idx_jjj);
            ome_num_step        = abs(metadata.scan_fin(idx_jjj) - metadata.scan_ini(idx_jjj))/data_ome_step_size;
            ome_num_step_calc   = CalcNumFramesGE(pfname_gefile);
            
            if (nframes_per_file ~= ome_num_step) || (nframes_per_file ~= ome_num_step_calc)
                is_good_ff_scan = false;
            end
        end
    end
    
    %%% NOW WRITE OUT PS FILES
    if is_good_ff_scan && isge
        if startsWith(froot_image_iii, froot_sample_name)
            ct  = ct + 1;
            
            fnum_first_layer    = metadata.det7_fnum(idx(1));
            nlayers_scan        = length(idx);
            
            disp(sprintf('GE_AD & FF-SCAN %s %d %d', ...
                froot_image_iii, fnum_first_layer, nlayers_scan));
            
            fname_ps_target{ct,1}   = sprintf('ps_midas_%s.txt', froot_image_iii);
            numlayers(ct,1)         = nlayers_scan;
            
            pfname_ps_template  = fullfile(pname_ps_template, fname_ps_template);
            pfname_ps_target    = fullfile(pname_ps_target, fname_ps_target{ct,1});
            
            cmdstr1 = sprintf('cp -v %s %s', pfname_ps_template, pfname_ps_target);
            status  = unix(cmdstr1);
            
            % str1    = sprintf('FileStem %s', gefile_prefix);
            % str2    = sprintf('StartFileNrFirstLayer %d', fnum_first_layer);
            
            cmdstr2 = sprintf('sed -i ''s/FileStem sam3_ff_155N/FileStem %s/g'' %s', ...
                froot_image_iii, pfname_ps_target);
            status  = unix(cmdstr2);
            
            cmdstr3 = sprintf('sed -i ''s/StartFileNrFirstLayer 269/StartFileNrFirstLayer %d/g'' %s', ...
                fnum_first_layer, pfname_ps_target);
            status  = unix(cmdstr3);
            
            cmdstr4 = sprintf('chmod -v 777 %s', ...
                pfname_ps_target);
            status  = unix(cmdstr4);
        end
    end
    pause(0.5)
end

pfname_midas_script = fullfile(pname_ps_target, fname_midas_script);
fid = fopen(pfname_midas_script, 'w');
for iii = 1:1:length(fname_ps_target)
    str1    = sprintf('/clhome/TOMO1/.MIDAS/MIDAS_V5_FarField_Layers %s 1 %d 1 11 orthrosnew parkjs@anl.gov\n', ...
        fname_ps_target{iii,1}, numlayers(iii,1));
    fprintf(fid, str1);
end
fclose(fid);
cmdstr  = sprintf('chmod -v 777 %s', ...
    pfname_midas_script);
status  = unix(cmdstr);