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

%%% Edits by Connor Horn on 7/20/18 to adapt to .cor32 format

iscor = isfield(ImagePars, 'iscor') && ImagePars.iscor;
iscor32 = isfield(ImagePars, 'iscor32') && ImagePars.iscor32;

if ~iscor && ~iscor32
    numimages   = length(ImagePars.fnumber);
    
    for i = 1:1:numimages
        % fname   = sprintf([ImagePars.fbase, '%05d.%s'], ImagePars.fnumber(i), ImagePars.fext);
        fname   = sprintf([ImagePars.fbase, '%0', num2str(ImagePars.numdigs), 'd.%s'], ImagePars.fnumber(i), ImagePars.fext);
        pfname{i,1} = fullfile(ImagePars.pname, fname);
    end
elseif iscor
    numimages   = length(ImagePars.fnumber);
    
    for i = 1:1:numimages
        % fname_base  = sprintf([ImagePars.fbase, '%05d.%s'], ImagePars.fnumber(i), ImagePars.fext);
        fname_base  = sprintf([ImagePars.fbase, '%0', num2str(ImagePars.numdigs), 'd.%s'], ImagePars.fnumber(i), ImagePars.fext);
        for j = 1:1:ImagePars.numframe
            %fname   = sprintf([fname_base, '.frame%d.cor'], j);
            fname   = sprintf([fname_base, '_frame_%d.cor'], j); % edited by CH on 8/3/17 to reflect format of batchcorr output by adding a '.'
            pfname{i,j} = fullfile(ImagePars.pname, fname);
        end
    end
    pfname  = pfname';
    pfname  = {pfname{:}}';
elseif iscor32
    numimages   = length(ImagePars.fnumber);
    
    for i = 1:1:numimages
        fname_base  = sprintf([ImagePars.fbase, '%0', num2str(ImagePars.numdigs), 'd.%s'], ImagePars.fnumber(i), ImagePars.fext);
        for j = 1:1:ImagePars.numframe
            %fname   = sprintf([fname_base, '.frame.%d.cor32'], j);
            fname   = sprintf([fname_base, '_frame_%d.cor32'], j);
            pfname{i,j} = fullfile(ImagePars.pname, fname);
        end
    end
    pfname  = pfname';
    pfname  = {pfname{:}}';
end

