% close all
clear all
clc

%%% look at GE image
fpath   = 'V:\Miller_July11\';
fname   = 'IN100_SP1_';
fnumber = 4611;

nn  = '00000';
nn(6-length(num2str(fnumber)):end) = num2str(fnumber);
bg  = NreadGE([fpath fname nn], 1);

img = 0*bg;
% 1397:1:1401
% 3500:1:3504
% 5503:1:5507
% 6868:1:6873
% 8850:1:8854
% 8876:1:8880
% 8904:1:8908
% 8914:1:8919
% 8968-9017
for i = 9003
    fnumber = i;
    fname   = 'IN100_SP1_';
    nn  = '00000';
    nn(6-length(num2str(fnumber)):end) = num2str(fnumber);
    im  = NreadGE([fpath fname nn], 1);
    im  = rot90(im-bg);
    img = img + im;
    
%     figure(i)
%     imagesc(im, [-20 3000]); colorbar
%     title(nn);
%     axis square
%     max(img(:))
    % pause()
end

figure(i)
imagesc(img, [-20 2000]); colorbar
title('sum');
max(img(:))
axis square