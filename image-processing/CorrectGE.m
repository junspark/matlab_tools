function [] = CorrectGE(pname_input, froot_img, fnum_img, froot_dark, fnum_dark, pname_output, keep_frames, ge_id, num_digit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprintf_cmd     = sprintf('%%s_%%0%dd.ge%%d', num_digit);
cmd_str         = sprintf('mkdir -v %s/%s', pname_output, froot_img);
[~, r]          = unix(cmd_str);
disp(r);

pout        = fullfile(pname_output, froot_img);

if ge_id == 1
    pfbad   = '/home/beams/S1IDUSER/mnt/s1a/misc/BatchCorrection/GE1Bad.img';
elseif ge_id == 2
    pfbad   = '/home/beams/S1IDUSER/mnt/s1a/misc/BatchCorrection/GE2Bad.img';
elseif ge_id == 3
    pfbad   = '/home/beams/S1IDUSER/mnt/s1a/misc/BatchCorrection/GE3Bad.img';
elseif ge_id == 4
    pfbad   = '/home/beams/S1IDUSER/mnt/s1a/misc/BatchCorrection/GE4Bad.img';
else
    disp('no bad file defined; check ge_id ...')
    return
end
im_bad  = NreadGE(pfbad, 1);

badpixel_type_2         = find(im_bad == 2);
badpixel_type_rem_one   = rem(im_bad, 2) == 1;

%%% GET DARK
fname_dark  = sprintf(sprintf_cmd, froot_dark, fnum_dark, ge_id);
pfname_dark = fullfile(pname_input, fname_dark);
num_frames = CalcNumFramesGE(pfname_dark);
im_dark = zeros(2048, 2048);
for i = 1:1:num_frames
    im_dark = im_dark + double(NreadGE(pfname_dark, i));
end
im_dark = im_dark./num_frames;

parpool(20);
parfor i = fnum_img
    fname_img   = sprintf(sprintf_cmd, froot_img, i, ge_id);
    pfname_img  = fullfile(pname_input, fname_img);
    
    fname_sum   = sprintf('%s.sum', fname_img);
    pfname_sum  = fullfile(pout, fname_sum);

    disp(sprintf('working on %s ...', fname_img));
    num_frames  = CalcNumFramesGE(pfname_img);
    
    im_imgi = zeros(2048, 2048);
    for j = 1:1:num_frames
        im_imgj = double(NreadGE(pfname_img, j));
        im_imgi = im_imgi + im_imgj;
        if keep_frames
            fname_cor   = sprintf('%s.frame.%d.cor', fname_img, j);
            pfname_cor  = fullfile(pout, fname_cor);
            
            im_imgj = im_imgj - im_dark;
            im_imgj = RemoveBadPixelGE(im_imgj, badpixel_type_2, badpixel_type_rem_one);
            
            disp(sprintf('saving corrected frame %d of %d to %s ...', j, num_frames, fname_cor));
            WriteSUM(pfname_cor, im_imgj);
        end
    end
    im_imgi = im_imgi - im_dark*num_frames;
    im_imgi = RemoveBadPixelGE(im_imgi, badpixel_type_2, badpixel_type_rem_one);
    disp(sprintf('saving corrected sum image to %s ...', fname_sum));
    WriteSUM(pfname_sum, im_imgi);
end
delete(gcp)