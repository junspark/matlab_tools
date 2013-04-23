function [xdata, ydata] = ReadEDDData(pfname, varargin)
% PlotSPF - plot point data on sphere
%
%   USAGE:
%
%   [xdata, ydata] = ReadEDDData(pfname)
%   [xdata, ydata] = ReadEDDData(pfname, varargin)
%
%   INPUT:
%
%   pfname
%       full file name of the EDD data file
%
%   These arguments can be followed by a list of
%   parameter/value pairs. Options are:
%
%   'Channel2Energy' line equation to convert channel number to energy
%   (keV)
%
%   OUTPUT:
%
%   xdata
%       x-coordinate of the data (in keV if conversion is provided)
%
%   ydata
%       y-coordinate of the data (counts)

% default options
optcell = {...
    'Channel2Energy', [1 1], ...
    };

% update option
opts    = OptArgs(optcell, varargin);

fid = fopen(pfname, 'r');
while ~feof(fid)
    lindata = fgetl(fid);
    if ~isempty(strfind(lindata, 'CHANNELS:'))
        numdata = str2num(lindata(20:end));
        ydata   = zeros(numdata,1);
        xdata   = 1:1:numdata;
    end
    
    if strcmp(lindata, 'DATA: ')
        for j = 1:1:numdata
            ydata(j,1)  = str2num(fgetl(fid));
        end
        break
    end
end
fclose(fid);

xdata   = opts.Channel2Energy(1)*xdata + opts.Channel2Energy(2);