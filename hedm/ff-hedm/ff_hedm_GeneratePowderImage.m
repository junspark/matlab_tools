function status = ff_hedm_GeneratePowderImage(pname, fstem, imnum, fext, fout)
% ff_hedm_GeneratePowderImage - generates synthetic powder image from hedm
% data. it is assumed that the input files are already summed in omega
% using batch correction routines.
%
% Input:
%
%   pname: path where ge sum files exist
%
%   fstem: stem of the ge sum files
%
%   imnum: n x 1 array of file numbers
%
%   fext: extension of the input files
%
%   fout: file name for the output file (synthetic powder patterns is
%   output in .sum format)
%

imsum   = zeros(2048,2048);
for i = 1:1:length(imnum)
    fname   = [fstem, '_', sprintf('%05d', imnum(i)), fext];
    pfname  = fullfile(pname, fname);
    
    imi = ReadSUM(pfname);
    
    imagesc(log(abs(imi)))
    axis equal
    pause(1)
    
    imsum   = imsum + imi;
end

pfname  = fullfile(pname, fout);
status  = WriteSUM(pfname, imsum);