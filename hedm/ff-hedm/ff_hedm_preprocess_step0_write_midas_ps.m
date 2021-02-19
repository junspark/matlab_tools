clear all
close all
clc

addpath(genpath('/home/beams0/PARKJS/matlab/matlab_tools'));

geprefix    = 'mishra_hr2_s1_load';
loadnum     = 0:1:14;
% geprefix	= 'mishra_hr1_s1_load';
% loadnum     = 0:1:10;
% geprefix	= 'mishra_hr3_s1_load';
% loadnum     = 0:1:15;

pfname_metadata     = '/home/beams/S1IDUSER/new_data/mpe1_oct20/fastpar_mpe1_oct20_FF.par';

pname_ps_template   = '/home/beams/S1IDUSER/mnt/orthros/mpe1_oct20_midas/ff/mishra_CuHEA_hcp/';
fname_ps_template   = 'ps_Cu_HEA_hcp_template.txt';
pname_ps_target     = pname_ps_template;
fname_midas_script  = 'run_midas.out';

metadata = readtable(pfname_metadata, 'FileType', 'text');


for iii = 1:1:length(loadnum)
    gefile_prefix   = sprintf('%s%d_ff', geprefix, loadnum(iii));
    imnum           = find(strcmp(metadata.Var22, gefile_prefix));
    
    fname_ps_target{iii,1}  = sprintf('%s.txt', gefile_prefix);
    numlayers(iii,1)        = length(imnum);
    
    pfname_ps_template  = fullfile(pname_ps_template, fname_ps_template);
    pfname_ps_target    = fullfile(pname_ps_target, fname_ps_target{iii,1});
    
    cmdstr1 = sprintf('cp -v %s %s', pfname_ps_template, pfname_ps_target);
    status  = unix(cmdstr1);
    
    str1    = sprintf('FileStem %s', gefile_prefix);
    str2    = sprintf('StartFileNrFirstLayer %d', imnum(1));
    
    cmdstr2 = sprintf('sed -i ''s/FileStem mishra_hr2_s1_load0_ff/FileStem %s/g'' %s', gefile_prefix, pfname_ps_target);
    status  = unix(cmdstr2);
    
    cmdstr3 = sprintf('sed -i ''s/StartFileNrFirstLayer 23/StartFileNrFirstLayer %d/g'' %s', imnum(1), pfname_ps_target);
    status  = unix(cmdstr3);
    
    cmdstr4 = sprintf('chmod -v 777 %s', pfname_ps_target);
    status  = unix(cmdstr4);
end

pfname_midas_script = fullfile(pname_ps_target, fname_midas_script);
fid = fopen(pfname_midas_script, 'w');
for iii = 1:1:length(fname_ps_target)
    str1    = sprintf('/clhome/TOMO1/.MIDAS/MIDAS_V5_FarField_Layers %s 1 %d 1 11 orthrosnew parkjs@anl.gov\n', ...
        fname_ps_target{iii,1}, numlayers(iii,1));
    fprintf(fid, str1);
end
fclose(fid);
% cmdstr2 = sprintf('chmod -v 777 %s', pfname_midas_script);
% status  = unix(cmdstr2);
