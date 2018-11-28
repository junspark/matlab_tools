function num_frames = CalcNumFramesGE(pfname)
% CalcNumFramesGE - Calculates the numbers of frames in a GE file with
% multiple frames.
%
%   USAGE:
%
%   num_frames = CalcNumFramesGE(pfname)
%
%   INPUT:
%
%   pfname          path where GE file is located.
%
%   OUTPUT:
%
%   num_frames      number of frames in the file.

if exist(pfname, 'file')
    num_frames  = dir(pfname);
    num_frames  = (num_frames.bytes-8192)/2048/2048/2;
else
    disp(sprintf('%s does not exist.', pfname))
    num_frames  = -1;
end