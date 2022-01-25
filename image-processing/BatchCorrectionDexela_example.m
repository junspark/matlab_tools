clear all
close all
clc

path_bkg    = '/home/beams/S1IDUSER/mnt/s1c/shastri_jul19/dexela';
path_image  = '/home/beams/S1IDUSER/mnt/s1c/shastri_jul19/dexela';

%%%% NEED TO CHANGE VARIABLE NAMES 
bkg_num     = 1206:1215;               %%% MAKE SURE EXPOSURE TIME IS SAME FOR THE BACKGROUND IMAGE AND DATA TO CORRECT. TYPICALLY dark_before TAKEN BEFORE EACH SHOULD WORK.
root_bkg    = 'dark_before';   %%% NOTE "_" AT THE END OF THE PREFIX
im_num      = [1216:1515, 1556:1855];
root_image  = 'AlCe_sam2_load0'; %%% NOTE "_" AT THE END OF THE PREFIX
path_output = '/home/s1b/__eval/projects_parkjs/shastri_jul19/dunand_AlCe_alloy';   %%% CHANGE OUTPUT DIRECTORY TO REFLECT YOUR COMPUTING ENV
%%%%

for i = 1:1:length(bkg_num)
    fname   = sprintf('%s_%06d.tif', root_bkg, bkg_num(i))
    pfname  = fullfile(path_bkg, root_image, fname);
    bkgi    = double(ReadDexela(pfname));
    if i == 1
        bkg     = bkgi;
    else
        bkg     = bkgi + bkg;
    end
end
bkg = bkg./length(bkg_num);
% imagesc(bkg)

img = bkg.*0;
parfor i = 1:1:length(im_num)
    fname   = sprintf('%s_%06d.tif', root_image, im_num(i))
    pfname  = fullfile(path_bkg, root_image, fname);
    imgi    = double(ReadDexela(pfname));
    
    img = imgi + img;
end
img = img - length(im_num)*bkg;
img = img./length(im_num);
img = img + abs(min(img(:)));
imagesc(img)

fname_out	= sprintf('%s_%06d_%06d.sum.tif', root_image, im_num(1), im_num(end));
pfname_out  = fullfile(path_output, fname_out);
% imwrite(img, pfname_out, 'Compression', 'none')
fid = fopen(pfname_out, 'w');
status  = fwrite(fid, img, 'int32');
status  = fclose(fid);