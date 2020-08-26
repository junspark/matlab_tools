clear all
close all
clc

addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

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

pname_lineout0  = 'C:\Users\parkjs\Documents\MATLAB\work\daymond_jun20_analysis\TEM_foil\lineout';

for iii = 1:1:length(fprefix)
    pname_lineout   = fullfile(pname_lineout0, fprefix{iii});
    for jjj = 1:1:length(fnum{iii})
        fname_MIDAS_csv     = sprintf('%s_%06d.raw_sum.csv', fprefix{iii}, fnum{iii}(jjj));
        pfname_MIDAS_csv    = fullfile(pname_lineout, fname_MIDAS_csv);
        
        fname_GSAS2_chi     = sprintf('%s_%06d.raw_sum.chi', fprefix{iii}, fnum{iii}(jjj))
        pfname_GSAS2_chi    = fullfile(pname_lineout, fname_GSAS2_chi);
        
        Convert_MIDAS_csv2GSAS2_chi(pfname_MIDAS_csv, pfname_GSAS2_chi, fname_GSAS2_chi, 0)
    end
end
% 