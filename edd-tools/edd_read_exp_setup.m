function [run_id, exp_id, pfname_edf, d_pk, d_UB, d_LB, froot, pname_mat, save_mat] = edd_read_exp_setup(exp_setup_file)
% edd_read_exp_setup - reads EDD experiment setup values from markdown.
%
%   USAGE:
%
%   [run_id, exp_id, pfname_edf] = edd_read_exp_setup(exp_setup_file)
%   [run_id, exp_id, pfname_edf, d_pk, d_UB, d_LB, froot, pname_mat, save_mat] = ...
%       edd_read_exp_setup(exp_setup_file)
%
%   INPUT:
%
%   exp_setup_file
%       Full path to exp_setup.md (or similar text file) containing
%       key:value entries.
%
%   OUTPUT:
%
%   run_id
%       Run identifier string parsed from key 'run_id'.
%
%   exp_id
%       Experiment identifier string parsed from key 'exp_id'.
%
%   pfname_edf
%       EDF file path string parsed from key 'pfname_edf'.
%
%   d_pk, d_UB, d_LB
%       Numeric row vectors parsed from optional keys 'd_pk', 'd_UB', and
%       'd_LB'. Empty if keys are not present.
%
%   froot, pname_mat
%       Strings parsed from optional keys 'froot' and 'pname_mat'. Empty if
%       keys are not present.
%
%   save_mat
%       Logical parsed from optional key 'save_mat'. Defaults to false if key
%       is not present.
%
%   REQUIRED FILE FORMAT:
%
%   run_id: <value>
%   exp_id: <value>
%   pfname_edf: <value>
%   d_pk: <space-separated numeric values>
%   d_UB: <space-separated numeric values>
%   d_LB: <space-separated numeric values>
%   froot: <value>
%   pname_mat: <value>
%   save_mat: <true|false|1|0|yes|no>

if ~isfile(exp_setup_file)
    error('Missing setup file: %s', exp_setup_file);
end

exp_setup_text = fileread(exp_setup_file);

% Defaults for optional values.
d_pk = [];
d_UB = [];
d_LB = [];
froot = '';
pname_mat = '';
save_mat = false;

% Read required key:value entries from exp_setup.md.
% Use [^\n]+ instead of .+ to exclude newlines from capture group.
run_id_token = regexp(exp_setup_text, '(?m)^\s*run_id\s*:\s*([^\n]+)\s*$', 'tokens', 'once');
exp_id_token = regexp(exp_setup_text, '(?m)^\s*exp_id\s*:\s*([^\n]+)\s*$', 'tokens', 'once');
pfname_token = regexp(exp_setup_text, '(?m)^\s*pfname_edf\s*:\s*([^\n]+)\s*$', 'tokens', 'once');

if isempty(run_id_token) || isempty(exp_id_token) || isempty(pfname_token)
    error('Invalid exp_setup.md format. Required keys: run_id, exp_id, pfname_edf');
end

run_id = strtrim(run_id_token{1});
exp_id = strtrim(exp_id_token{1});
pfname_edf = strtrim(pfname_token{1});

% Read optional key:value entries when present.
d_pk_str = read_optional_key(exp_setup_text, 'd_pk');
if ~isempty(d_pk_str)
    d_pk = sscanf(d_pk_str, '%f').';
end

d_UB_str = read_optional_key(exp_setup_text, 'd_UB');
if ~isempty(d_UB_str)
    d_UB = sscanf(d_UB_str, '%f').';
end

d_LB_str = read_optional_key(exp_setup_text, 'd_LB');
if ~isempty(d_LB_str)
    d_LB = sscanf(d_LB_str, '%f').';
end

froot_str = read_optional_key(exp_setup_text, 'froot');
if ~isempty(froot_str)
    froot = strtrim(froot_str);
end

pname_mat_str = read_optional_key(exp_setup_text, 'pname_mat');
if ~isempty(pname_mat_str)
    pname_mat = strtrim(pname_mat_str);
end

save_mat_str = read_optional_key(exp_setup_text, 'save_mat');
if ~isempty(save_mat_str)
    save_mat = parse_logical_value(save_mat_str, 'save_mat');
end
end

function value = read_optional_key(text_blob, key_name)
pattern = sprintf('(?m)^\\s*%s\\s*:\\s*([^\\n]+)\\s*$', regexptranslate('escape', key_name));
token = regexp(text_blob, pattern, 'tokens', 'once');
if isempty(token)
    value = '';
else
    value = strtrim(token{1});
end
end

function value = parse_logical_value(raw_value, key_name)
normalized = lower(strtrim(raw_value));
if any(strcmp(normalized, {'true', '1', 'yes', 'y'}))
    value = true;
elseif any(strcmp(normalized, {'false', '0', 'no', 'n'}))
    value = false;
else
    error('Invalid logical value for key ''%s'': %s', key_name, raw_value);
end
end
