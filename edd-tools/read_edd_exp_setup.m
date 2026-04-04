function [run_id, exp_id, pfname_edf] = read_edd_exp_setup(exp_setup_file)
% Parse experiment setup values from a simple key:value markdown file.

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
