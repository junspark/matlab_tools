clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ge_id           = 1:4;
num_digit       = 6;
keep_frames     = 0;

pname_input     = '/home/beams/S1IDUSER/mnt/s1c/xu_nov16/ge';
pname_output    = '/home/beams/S1IDUSER/mnt/orthros/xu_nov16_bc/';

froot_dark      = 'dark_before';
fnum_dark       = 3099;

froot_img       = 'HTUPS_HT2_400C_go';
fnum_img        = 2540:3098;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:length(ge_id)
    CorrectGE(pname_input, froot_img, fnum_img, froot_dark, fnum_dark, pname_output, keep_frames, ge_id(i), num_digit)
end