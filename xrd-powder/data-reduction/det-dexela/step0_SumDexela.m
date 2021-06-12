clear all
close all
clc

addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

% 21 - fname
% 22 - scan begin 
% 23 - scan end
proot_out   = '/home/beams/S1IDUSER/mnt/orthros/internal_may20_bc';

proot   = '/home/beams/S1IDUSER/mnt/orthros/internal_may20_data/internal_may20/dexela';
froot   = 'joha_texture';
fdark       = 'dark_before';
fnum_dark   = 188:208;

% proot   = '/home/beams/S1IDUSER/mnt/orthros/internal_may20_data/internal_may20/dexela';
% froot   = 'CeO2_0pt5s';
% fdark       = 'dark_0pt5s';
% fnum_dark   = 146:166;

pmeta   = '/home/beams/S1IDUSER/new_data/internal_may20';
fmeta   = 'fastpar_internal_may20_Tomo.par';
pfmeta  = fullfile(pmeta, fmeta);

opts        = detectImportOptions(pfmeta,'FileType','text');
metadata    = readtable(pfmeta, opts);

idx = find(ismember(metadata.Var21, froot));

for iii = 1:1:length(fnum_dark)
    fname_img   = sprintf('%s_%06d.tif', fdark, fnum_dark(iii));
    pfname_img  = fullfile(proot, fname_img);
    
    imiii  = double(imread(pfname_img));
    
    if iii == 1
        dark    = imiii;
    else
        dark   = dark + imiii;
    end
end
dark    = dark./length(fnum_dark);

for iii = 1:1:length(idx)
    metadata.Var21(idx(iii))
    if metadata.Var23(idx(iii)) >= metadata.Var22(idx(iii))
        imnum   = metadata.Var22(idx(iii)):metadata.Var23(idx(iii));
        for jjj = 1:1:length(imnum)
            fname_img   = sprintf('%s_%06d.tif', froot, imnum(jjj));
            pfname_img  = fullfile(proot, fname_img);
            
            imjjj   = double(imread(pfname_img));
            imjjj2  = Tiff(pfname_img, 'r');
%             figure, 
%             imagesc(imjjj2.read)
%             figure, 
%             imagesc(imjjj)
%             return
            
            if jjj == 1
                imout   = imjjj;
            else
                imout   = imout + imjjj;
            end
        end
        imout   = imout./length(imnum);
        imout   = imout - dark;
        
        imout2  = uint16(imout);
        
        fout    = sprintf('%s_sum_%06d_%06d.sum.tif', froot, imnum(1), imnum(end));
        pfout   = fullfile(proot_out, fout);
        
        WriteDexela(pfout, imout2);
    end
end