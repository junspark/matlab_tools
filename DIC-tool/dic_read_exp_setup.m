function [pname, froot, fext, ndigits, fini, ffin, finc, ...
    delta_h, delta_v, Hctr, Vctr, ws_h, ws_v, pix2mm, gauge_length] = dic_read_exp_setup(exp_info_file)
% dic_read_exp_setup - reads DIC experiment setup values from markdown.
%
%   USAGE:
%
%   [pname, froot, fext, ndigits, fini, ffin, finc, ...
%       delta_h, delta_v, Hctr, Vctr, ws_h, ws_v, pix2mm, gauge_length] = dic_read_exp_setup(exp_info_file)
%
%   INPUT:
%
%   exp_info_file
%       Full path to dic_exp_info.md (or similar text file) containing
%       key:value entries.
%
%   OUTPUT:
%
%   pname       Path to directory containing DIC images.
%   froot       File root name (prefix before the frame number).
%   fext        File extension (e.g. tif).
%   ndigits     Number of digits used in the frame number.
%   fini        First frame number.
%   ffin        Last frame number.
%   finc        Frame number increment.
%   delta_h     Horizontal spacing between control points (pixels).
%   delta_v     Vertical spacing between control points (pixels).
%   Hctr        Horizontal center of ROI (pixels).
%   Vctr        Vertical center of ROI (pixels).
%   ws_h        Half-width of ROI in horizontal direction (control points).
%   ws_v        Half-width of ROI in vertical direction (control points).
%   pix2mm        Scale factor (mm per pixel).
%   gauge_length  Gauge length (mm).
%
%   REQUIRED FILE FORMAT:
%
%   pname: <value>
%   froot: <value>
%   fext: <value>
%   ndigits: <value>
%   fini: <value>
%   ffin: <value>
%   finc: <value>
%   delta_h: <value>
%   delta_v: <value>
%   Hctr: <value>
%   Vctr: <value>
%   ws_h: <value>
%   ws_v: <value>
%   pix2mm: <value>
%   gauge_length: <value>    (mm)

if ~isfile(exp_info_file)
    error('Missing setup file: %s', exp_info_file);
end

txt = fileread(exp_info_file);

% Helper to parse a string value
    function val = str_val(key)
        tok = regexp(txt, sprintf('(?m)^\\s*%s\\s*:\\s*([^\\n]+)\\s*$', key), 'tokens', 'once');
        if isempty(tok)
            error('Key "%s" not found in %s', key, exp_info_file);
        end
        val = strtrim(tok{1});
    end

% Helper to parse a numeric value
    function val = num_val(key)
        val = str2double(str_val(key));
    end

pname   = str_val('pname');
froot   = str_val('froot');
fext    = str_val('fext');
ndigits = num_val('ndigits');
fini    = num_val('fini');
ffin    = num_val('ffin');
finc    = num_val('finc');
delta_h = num_val('delta_h');
delta_v = num_val('delta_v');
Hctr    = num_val('Hctr');
Vctr    = num_val('Vctr');
ws_h    = num_val('ws_h');
ws_v         = num_val('ws_v');
pix2mm       = num_val('pix2mm');
gauge_length = num_val('gauge_length');

end
