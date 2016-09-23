function pfname = GenerateGEpfname(ImagePars)
% GenerateGEpfname - generate file names for GE image files
%
%   USAGE:
%
%   pfname = GenerateGEpfname(ImagePars)
%
%   INPUT:
%
%   ImagePars
%       Structure array that contains information about the file name
%       patterns. Must contain pname (path name), fbase (file name), fnumber(series of file
%       number), fext (file extension, i.e ge1)
%
%   OUTPUT:
% 
%   pfname
%       cell structure with the full file names
% 

if ~isfield(ImagePars, 'iscor') || ~ImagePars.iscor
    numimages   = length(ImagePars.fnumber);
    
    for i = 1:1:numimages
        % fname   = sprintf([ImagePars.fbase, '%05d.%s'], ImagePars.fnumber(i), ImagePars.fext);
        fname   = sprintf([ImagePars.fbase, '%0', num2str(ImagePars.numdigs), 'd.%s'], ImagePars.fnumber(i), ImagePars.fext);
        pfname{i,1} = fullfile(ImagePars.pname, fname);
    end
else
    numimages   = length(ImagePars.fnumber);
    
    for i = 1:1:numimages
        % fname_base  = sprintf([ImagePars.fbase, '%05d.%s'], ImagePars.fnumber(i), ImagePars.fext);
        fname_base  = sprintf([ImagePars.fbase, '%0', num2str(ImagePars.numdigs), 'd.%s'], ImagePars.fnumber(i), ImagePars.fext);
        for j = 1:1:ImagePars.numframe
            fname   = sprintf([fname_base, '.frame%d.cor'], j);
            pfname{i,j} = fullfile(ImagePars.pname, fname);
        end
    end
    pfname  = pfname';
    pfname  = {pfname{:}}';
end

