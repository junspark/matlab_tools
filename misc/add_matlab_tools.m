function status = add_matlab_tools()
status  = 0;
if isunix
    addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
    status  = 1;
elseif ispc
    pcname = getenv('COMPUTERNAME'); % for windows
    if strcmpi(pcname, 'sec1parks')
        addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));
        status  = 1;
    elseif strcmpi(pcname, 'riesling')
        addpath(genpath('W:\__eval\projects_parkjs\matlab_tools'));
        status  = 1;
    else
        disp(sprintf('%s is an unknown host', pcname))
    end
end