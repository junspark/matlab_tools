clear all
close all
clc

pfname_in   = 'W:\__eval\li_march13\HTUPS_Irr_00201.ge3';
pfname_out  = 'W:\__eval\li_march13\HTUPS_Irr_NoEmptyFrame00201.ge3';

nFrames = 182;
nFrames_Empty   = 2;

%%% WRITING PART
for i = 1:1:(nFrames-nFrames_Empty)
    img(:,:,i)  = NreadGE(pfname_in, (i+nFrames_Empty) );
    
    % figure(1)
    % imagesc(log(abs(double(img))))
    % caxis([0 log(500)])
    % colorbar vert
    % axis square tight
    % pause
end
status  = NwriteGE(pfname_out, img)

%%% CHECKING PART
% for i = 1:1:(nFrames-nFrames_Empty)
%     img_in	= NreadGE(pfname_in, (i+nFrames_Empty) );
%     img_out = NreadGE(pfname_out, i);
%     
%     figure(1)
%     subplot(1,2,1)
%     imagesc(log(abs(double(img_in))))
%     caxis([0 log(500)])
%     colorbar vert
%     axis square tight
%     
%     subplot(1,2,2)
%     imagesc(log(abs(double(img_out))))
%     caxis([0 log(500)])
%     colorbar vert
%     axis square tight
%     
%     diff(i) = sum(img_in(:) - img_out(:));
% end
% 
% figure,
% plot(diff)