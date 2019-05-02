function output_vals = ConvertMIDASToHEXRD(input_vals, varargin)
% NEED TO CHECK THE ETA CONVERSIONS 
% NEED TO IMPLEMENT OTHER CONVERSIONS
% ConvertMIDASToHEXRD Converts various parameters FROM MIDAS TO HEXRD
% CONVENTION
%
%   INPUT:
%
%   input_vals
%       values to be converted
%
%   These arguments can be followed by a list of
%   parameter/value pairs. Options are:
%
%   'ObjectToConvert'   eta only (eta is default)
%   'Units'             deg or rad (deg is default)
%
%   OUTPUT:
%   
%   output_vals
%       converted values. if invalid input_vals then returns nan

% default options
optcell = {...
    'ObjectToConvert', 'eta', ...
    'Units', 'deg', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

output_vals = nan;
if strcmpi(opts.ObjectToConvert, 'eta')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% HEXRD ETA CONVERSION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(opts.Units, 'rad')
        input_vals  = rad2deg(input_vals);
    elseif ~strcmpi(opts.Units, 'deg')
        disp('angular unit conversion not implemented yet!')
        return
    end
    
    output_vals = -input_vals+90;
    idx         = output_vals > 180;
    if sum(idx) > 0
        output_vals(idx)    = -input_vals(idx)-270;
    end
else
    disp('object conversion not implemented yet!')
    return
end