function pardata = ReadSpecPerFrameParFile(pfname)
% ReadSpecPerFrameParFile - read per frame par file generated from spec for
% WAXS data
%
%   INPUT:
%
%   fname
%       name of the spec generated per frame par file name
%
%   OUTPUT:
%
%   pardata
%       spec par file data in struct arrary.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The input arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'Version'               version of spec metadata file (par file).
%                           example: (coratella_feb15). default value is
%                           'none'. 'mpe_standard' is the new MPE metadata
%                           format for 1-BM-B EDD / 6-BM-A EDD / 1-ID data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all
% clc
% 
% pname   = '/home/beams/S1IDUSER/new_data/fmeng_apr18/';
% fname   = 'per_frame_waxs_hydra_exposures2.par';
% 
% pfname  = fullfile(pname, fname);

% fmtstring0  = '%s %s %d %s %d %f %f %f %s %d %d %f %f';

% for i = 1:1:numchs
%             fmtstring   = [fmtstring, ' %f'];
% end
%  p date(),sprintf("%15.8f", epoch_time),elapsed_time,Iring, \
%             GE_fname,GE_fnum3,GE_Nframe,GE_tframe,GE_tframe,\
%             timelist, temp1list, temp2list, temp3list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(pfname, 'r');
ct  = 0;
idx_start_of_per_frame_data = 14;
while ~feof(fid)
    ct  = ct + 1;
    lindata = fgetl(fid);
    lindata = strsplit(lindata);
    
    day{ct}     = lindata{1};
    month{ct}   = lindata{2};
    date(ct)    = str2num(lindata{3});
    time{ct}    = lindata{4};
    year{ct}    = str2num(lindata{5});
    
    epoch_time{ct}  = str2num(lindata{6});
    integ_time{ct}  = str2num(lindata{7});
    Iring{ct}       = str2num(lindata{8});
    
    det_fname{ct}           = lindata{9};
    det_fnum(ct)            = str2num(lindata{10});
    det_frames_per_file(ct) = str2num(lindata{11});
    det_time_per_frame(ct)  = str2num(lindata{12});
    
    %%% GET PER FRAME DATA
    nframes = det_frames_per_file(ct);
    
    %%% SET 0
    idx1    = idx_start_of_per_frame_data;
    idx2    = idx_start_of_per_frame_data + nframes - 1;
    
    per_frame_data      = str2mat(lindata{idx1:idx2});
    per_frame_data0{ct} = str2num(per_frame_data);
    
    %%% SET 1
    idx1    = idx2 + 1;
    idx2    = idx1 + nframes - 1;
    
    per_frame_data      = str2mat(lindata{idx1:idx2});
    per_frame_data1{ct} = str2num(per_frame_data);
    
    %%% SET 2
    idx1    = idx2 + 1;
    idx2    = idx1 + nframes - 1;
    
    per_frame_data      = str2mat(lindata{idx1:idx2});
    per_frame_data2{ct} = str2num(per_frame_data);
    
    %%% SET 3
    idx1    = idx2 + 1;
    idx2    = idx1 + nframes - 1;
    
    per_frame_data      = str2mat(lindata{idx1:idx2});
    per_frame_data3{ct} = str2num(per_frame_data);
end
fclose(fid);

pardata.day     = day;
pardata.month	= month;
pardata.date	= date;
pardata.time	= time;
pardata.year	= year;

pardata.epoch_time	= epoch_time;
pardata.integ_time	= integ_time;
pardata.Iring       = Iring;

pardata.det_fname           = det_fname;
pardata.det_fnum            = det_fnum;
pardata.det_frames_per_file = det_frames_per_file;
pardata.det_time_per_frame	= det_time_per_frame;

pardata.per_frame_data0 = per_frame_data0;
pardata.per_frame_data1 = per_frame_data1;
pardata.per_frame_data2 = per_frame_data2;
pardata.per_frame_data3 = per_frame_data3;