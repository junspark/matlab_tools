function [run_id, exp_id, pfname_edf] = edd_read_exp_setup(exp_setup_file)
% edd_read_exp_setup - reads EDD experiment setup values from markdown.
%
%   USAGE:
%
%   [run_id, exp_id, pfname_edf] = edd_read_exp_setup(exp_setup_file)
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
%   REQUIRED FILE FORMAT:
%
%   run_id: <value>
%   exp_id: <value>
%   pfname_edf: <value>

if ~isfile(exp_setup_file)
    error('Missing setup file: %s', exp_setup_file);
end

exp_setup_text = fileread(exp_setup_file);

% Read required key:value entries from exp_setup.md.
run_id_token = regexp(exp_setup_text, '(?m)^\s*run_id\s*:\s*(.+)\s*$', 'tokens', 'once');
exp_id_token = regexp(exp_setup_text, '(?m)^\s*exp_id\s*:\s*(.+)\s*$', 'tokens', 'once');
pfname_token = regexp(exp_setup_text, '(?m)^\s*pfname_edf\s*:\s*(.+)\s*$', 'tokens', 'once');

if isempty(run_id_token) || isempty(exp_id_token) || isempty(pfname_token)
    error('Invalid exp_setup.md format. Required keys: run_id, exp_id, pfname_edf');
end

run_id = strtrim(run_id_token{1});
exp_id = strtrim(exp_id_token{1});
pfname_edf = strtrim(pfname_token{1});
end
