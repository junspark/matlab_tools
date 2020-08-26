clear all
close all
clc

addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

% neta            = 72;
% nrbins          = 1055;
% Dsam            = 1476035.903248703340;

pname_pattern   = '/home/beams/S1IDUSER/mnt/orthros/daymond_jun20_data/daymond_jun20/pilatus';
pname_lineout   = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/daymond_jun20_analysis/TEM_foil/lineout';

fprefix{1}  = 'samTopLeft';
fnum{1}     = 23:38;

fprefix{2}  = 'samTopMid';
fnum{2}     = 41:56;

fprefix{3}  = 'samTopRight';
fnum{3}     = 59:74;

fprefix{4}  = 'samMidLeft';
fnum{4}     = 77:92;

fprefix{5}  = 'samMidMid';
fnum{5}     = 95:110;

fprefix{6}  = 'samMidRight';
fnum{6}     = 113:128;

fprefix{7}  = 'samBotLeft';
fnum{7}     = 131:146;

fprefix{8}  = 'samBotMid';
fnum{8}     = 149:164;

fprefix{9}  = 'samBotRight';
fnum{9}     = 167:182;

delete(gcp)
parpool(4)
for ii = 1:1:length(fprefix)
    for iii = 1:1:length(fnum{ii})
        fname_pattern   = sprintf('%s_%06d.raw', fprefix{ii}, fnum{ii}(iii));
        pfname_pattern  = fullfile(pname_pattern, fname_pattern);
        
        cmd_str1    = sprintf('!  ~/opt/MIDAS/FF_HEDM/bin/Integrator ps_integrate.txt %s', pfname_pattern);
        eval(cmd_str1)
        
        pname_target = fullfile(pname_lineout, fprefix{ii});
        if ~isdir(pname_target)
            mkdir(pname_target)
        end
        
        pause(1)
        cmd_str2    = sprintf('! mv -v %s/%s_* %s/%s', pname_lineout, fname_pattern, pname_lineout, fprefix{ii})
        eval(cmd_str2)
    end
end
delete(gcp)
