clear all
close all
clc

addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

% def daymond_blister_radiograph '{
%     local _iX _nX _dX _Xo _Xi _nowX
%     local _iY _nY _dY _Yo _Yi _nowY
%     local _et
%
%     switch_to_TomoDet
%     switch_to_TomoBeam
%
%     _dX=2;
%     _dY=0.9;
%     _nX=3;
%     _nY=4;
%     _Xi=-2.6;
%     _Yi=5.6;
%     for (_iX=0;_iX<_nX;_iX++) {
%         _nowX = _Xi + _iX*_dX;
%         for (_iY=0;_iY<_nY;_iY++) {
%             _nowY = _Yi + _iY*_dY;
%             umv samXE _nowX samYE _nowY
%
%             p _nowX, " ", _nowY
%             takeRadio "$0"
%             sleep(1)
%         }
%     }
% }'

pname_wf    = 'C:\Users\parkjs\Documents\MATLAB\work\daymond_jun20_data\tomo\broderick_dec16_S9_tomo_still_20x_L5';
fname_wf    = 'broderick_dec16_S9_tomo_still_20x_L5_036314.tif';
pfname_wf   = fullfile(pname_wf, fname_wf);
imdata_wf   = double(imread(pfname_wf));

pname   = 'C:\Users\parkjs\Documents\MATLAB\work\daymond_jun20_data\tomo';
froot   = 'daymond_blister_radiograph_0deg';
fnumi   = 47222;
fnum_list   = [ 
    47222
    47223
    47224
    47225
    47226
    47227
    47228
    47230
    47231
    47232
    47234
    47235
    ];

dX      = 2; dY  = 0.9;
nX      = 3; nY  = 4;
Xi      = -2.6; Yi  = 5.6;
pix2mm  = 1.172/1000;

% pTL = [175 25];
% pTR = [185 1885];
% pBL = [1020 25];
% pBR = [1020 1885];
fov_size_v  = (1020-185)*pix2mm;
fov_size_h  = (1885-25)*pix2mm;

overlap_v   = fov_size_v - dY;
overlap_h   = fov_size_h - dX;

overlap_v_pix   = round(overlap_v/pix2mm);
overlap_h_pix   = round(overlap_h/pix2mm);

numimg  = nX*nY;
ct  = 0;
for iii = 1:1:nX
    for jjj = 1:1:nY
        fnum    = fnum_list(ct+1);
        fname   = sprintf('%s_%06d.tif', froot, fnum);
        pfname  = fullfile(pname, fname);
        
        imdata  = imread(pfname);
        
        %%% NORMALIZE
        imdata  = double(imdata)./imdata_wf;
        
        %%% CROP IMAGE TO SLITS
        imdata  = imdata(185:1020, 25:1885);
        
        if jjj == 1
            imdata_merged_jjj   = imdata;
        else
            imdata_merged_jjj   = [imdata_merged_jjj; imdata(overlap_v_pix:end,:)];
        end
        ct      = ct + 1;
%         
%         figure(ct)
%         imagesc(imdata_merged_jjj)
%         colormap gray
%         caxis([0 1])
    end
    
    if iii == 1
        imdata_merged   = imdata_merged_jjj;
    else
        imdata_merged   = [imdata_merged imdata_merged_jjj(:,overlap_h_pix:end)];
    end
%     figure,
%     imagesc(imdata_merged)
%     return
end

figure,
imagesc(imdata_merged)
axis equal tight
colormap gray
caxis([0 0.7])
xticklabels(round(xticks*pix2mm*10)/10-3.6)
yticklabels(round(yticks*pix2mm*10)/10+5.1)
% xticklabels(round(xticks*pix2mm*10)/10)
% yticklabels(round(yticks*pix2mm*10)/10)
xlabel('X (mm)')
ylabel('Y (mm)')
