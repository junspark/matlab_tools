function [imdata, h] = ReadSMV(filename, varargin)
% [imdata, h] = ReadSMV(filename, vargin)
%   reads .smv image frames from Dexela

% clear all
% close all
% clc
% 
% frameno     = 1;
% filename    = './1512-11620-Dark-High-5x1000.0ms-Binx11.smv';

if nargin > 2
    error('Usage [imdata, h] = ReadSMV(filename, frame_number)');
elseif nargin == 2 && ~isnumeric(varargin{1})
    error('frame_number is a positive integer');
elseif varargin{1} < 1
    error('frame_number is a positive integer');
end

fp  = fopen(filename, 'r', 'n');

headerString    = fread(fp, 512, 'uint8=>char');
headerParts     = regexp(headerString(:,1)',';','split');

tmpParts        = regexp(headerParts{1},'=','split');
h.header_bytes  = str2double(tmpParts{end});

tmpParts    = regexp(headerParts{2},'=','split');
h.filename  = tmpParts{end};

h.filetype  = 'smv';

%h.dim: Dimension of the image. '1' is for a single image '2' for a stack of images
tmpParts    = regexp(headerParts{4},'=','split');
h.dim       = str2num(tmpParts{end});

%h.size1: The width of the image.
tmpParts    = regexp(headerParts{7},'=','split');
h.width     = str2num(tmpParts{end});

%h.size2: The height of the image. 
tmpParts    = regexp(headerParts{8},'=','split');
h.height    = str2num(tmpParts{end});

%h.size3: The number of images saved.
tmpParts    = regexp(headerParts{9},'=','split');
h.nframes   = str2num(tmpParts{end});

tmpParts    = regexp(headerParts{10},'=','split');
h.size4     = str2num(tmpParts{end}); % NOT USED

%read image data
if nargin == 1
    imdata = zeros(h.height, h.width, h.nframes);
    for i = 1:1:h.nframes
        img = fread(fp, [h.width, h.height], 'ushort');
        img = rot90(img);
        imdata(:,:,i)   = img;
    end
elseif nargin == 2
    if varargin{1} > h.nframes
        error('frame_number needs to be smaller than the number of frames in the image stack');    
        fclose(fp);
    else
        frameno = varargin{1};
        fseek(fp, h.header_bytes + (frameno-1)*2*h.height*h.width, 'bof');
        
        imdata = fread(fp, [h.width, h.height], 'ushort');
        imdata = rot90(imdata);
    end
end
fclose(fp);

end