function pardata = ReadPythonParFile(fname, varargin)
% ReadPythonParFile - read par file generated from python
%
%   INPUT:
%
%   fname
%       name of the python generated par file name
%
%   pypar_version
%       versions of the pypar file
%       'current' is the current state
%       
%
%   OUTPUT:
%
%   pardata
%       strip chart data in struct arrary.

% default options
optcell = {...
    'pypar_version', 'feb15', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

if strcmpi(opts.pypar_version, 'feb15')
    fmtstring       = '%s %s %d %s %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s %s %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f';
elseif strcmpi(opts.pypar_version, 'dec14')
    fmtstring       = '%s %s %d %s %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s %s %f, %f, %f, %f, %f';
end

fid     = fopen(fname, 'r','n');
header  = fgetl(fid);
textdata  = textscan(fid, fmtstring);
fclose(fid);

if strcmpi(opts.pypar_version, 'feb15')
    pardata.hdr                 = header;
    pardata.day                 = textdata{1};
    pardata.month               = textdata{2};
    pardata.date                = textdata{3};
    pardata.time                = textdata{4};
    pardata.year                = textdata{5};
    pardata.Iring               = textdata{6};
    pardata.ElapsedTime         = textdata{7};
    pardata.ic0                 = textdata{8};
    pardata.ic1                 = textdata{9};
    pardata.NULL                = textdata{10};
    pardata.HorzInstDeadTime    = textdata{11};
    pardata.HorzAveDeadTime     = textdata{12};
    pardata.HorzRealElapTime    = textdata{13};
    pardata.HorzLiveElapTime    = textdata{14};
    pardata.VertInstDeadTime    = textdata{15};
    pardata.VertAveDeadTime     = textdata{16};
    pardata.VertRealElapTime    = textdata{17};
    pardata.VertLiveElapTime    = textdata{18};
    pardata.IterationNumber     = textdata{19};
    pardata.PositionNumber      = textdata{20};
    pardata.FileIDNumber        = textdata{21};
    pardata.HorzFileName        = textdata{22};
    pardata.VertFileName        = textdata{23};
    pardata.samX                = textdata{24};
    pardata.samY                = textdata{25};
    pardata.samY2               = textdata{26};
    pardata.samZ                = textdata{27};
    pardata.HorzSlitB           = textdata{28};
    pardata.HorzSlitT           = textdata{29};
    pardata.HorzSlitI           = textdata{30};
    pardata.HorzSlitO           = textdata{31};
    pardata.HorzTh              = textdata{32};
    pardata.VertSlitB           = textdata{33};
    pardata.VertSlitT           = textdata{34};
    pardata.VertSlitI           = textdata{35};
    pardata.VertSlitO           = textdata{36};
    pardata.VertTh              = textdata{37};
elseif strcmpi(opts.pypar_version, 'dec14')
    pardata.hdr                 = header;
    pardata.day                 = textdata{1};
    pardata.month               = textdata{2};
    pardata.date                = textdata{3};
    pardata.time                = textdata{4};
    pardata.year                = textdata{5};
    pardata.Iring               = textdata{6};
    pardata.ElapsedTime         = textdata{7};
    pardata.ic0                 = textdata{8};
    pardata.ic1                 = textdata{9};
    pardata.NULL                = textdata{10};
    pardata.HorzInstDeadTime    = textdata{11};
    pardata.HorzAveDeadTime     = textdata{12};
    pardata.HorzRealElapTime    = textdata{13};
    pardata.HorzLiveElapTime    = textdata{14};
    pardata.samX                = textdata{15};
    pardata.samY                = textdata{16};
    pardata.samY2               = textdata{17};
    pardata.samZ                = textdata{18};
    pardata.HorzSlitB           = textdata{19};
    pardata.HorzSlitT           = textdata{20};
    pardata.HorzSlitI           = textdata{21};
    pardata.HorzSlitO           = textdata{22};
    pardata.HorzTh              = textdata{23};
    pardata.VertSlitB           = textdata{24};
    pardata.VertSlitT           = textdata{25};
    pardata.VertSlitI           = textdata{26};
    pardata.VertSlitO           = textdata{27};
    pardata.VertTh              = textdata{28};
end
