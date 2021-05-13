function imctrl = GSAS2_Read_imctrl(pfname)
% GSAS2_Read_imctrl - GSAS2 imctrl file reader
%
%   INPUT:
%
%   pfname
%       full pathname of the GSAS2 imctrl file.
%
%   OUTPUT:
%
%   imctrl_data
%       string array data of the the GSAS2 imctrl file.

% imctrl_data = readtable(pfname, 'FileType', 'text', 'ReadVariableNames', false, ...
%     'delimiter', ':' ,'ConsecutiveDelimitersRule' ,'join', ...
%     'NumHeaderLines', 0, ...
%     'TextType', 'char', ...
%     'ReadRowNames', true);

% options = detectImportOptions(pfname, 'FileType', 'text', 'Delimiter', ':')
% imctrl_data = readtable(pfname, options);

fid = fopen(pfname ,'r');
imctrl_data = textscan(fid, '%s', 'CollectOutput', true, 'Delimiter', ':', 'TextType', 'string');
fclose(fid);

imctrl_data = string(imctrl_data{:});

idx                 = find(contains(imctrl_data, 'LRazimuth'));
imctrl.LRazimuth    = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'azmthOff'));
imctrl.azmthOff     = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'outAzimuths'));
imctrl.outAzimuths  = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'IOtth'));
imctrl.IOtth        = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'outChannels'));
imctrl.outChannels  = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'center'));
imctrl.center       = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'wavelength'));
imctrl.wavelength   = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'distance'));
imctrl.distance     = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'tilt'));
imctrl.tilt         = eval(imctrl_data{idx+1});

idx                 = find(contains(imctrl_data, 'rotation'));
imctrl.rotation     = eval(imctrl_data{idx+1});

%%% SOME DERIVED QUANTITIES FROM IMCTRL
dazm                = (imctrl.LRazimuth(2) - imctrl.LRazimuth(1))/imctrl.outAzimuths;
azm_grid            = linspace(imctrl.LRazimuth(1)+dazm/2, imctrl.LRazimuth(2)-dazm/2, imctrl.outAzimuths);
idx_azm             = azm_grid >= 360;
azm_grid(idx_azm)   = azm_grid(idx_azm) - 360;

imctrl.derived_dazm     = dazm;
imctrl.derived_azm_grid = azm_grid;