% EDD data file reader (6BM setup)
%  read data into a single Matlab struct
%  supported configuration: Nov.2016, Mar.2017, Aug.2017
%
%  Rev.1.6 (2017/08/09)
%  + combine detector data into cell
%
%  Rev.1.5 (2017/03/23)
%  + bug fix
%
%  Rev.1.4 (2017/03/15)
%  + read motor_cam2 data
%  + read all motor position data
%
%  Rev.1.3 (2017/01/15)
%  + read 2nd detector data
%  + read folder
%
%  Rev.1.2 (2016/11/30)
%  + add keyence to monitor yr (yr motor loss step issues)
%
%  Rev.1.0 (2015/12/18)
%  + to be documented
%  + read multiple scans from single file (append mode)
%
% Copyright 2015-2017 Andrew Chuang (chuang.cp@gmail.com)
% $Revision: 1.6 $  $Date: 2017/08/09 $


function [edd]=readedd3_6bm(logtoopen)

if nargin ~= 1
    fprintf('\nUsage : [log]=readedd_6bm("FILEtoOPEN")\n');
    fprintf('\n');
    return;
end

% open and read the file into memory
if isdir(logtoopen)
    [fpath, fname, ~] = fileparts(logtoopen);
    fullname= fullfile(logtoopen,sprintf('%s.xy',fname));
    fid=fopen(fullname);if fid == -1, error('Can''t find/open the input file.'); end;
    step_counter = zeros(5000,1);
else
    fid=fopen(logtoopen);if fid == -1, error('Can''t open the input file.'); end;
    step_counter = zeros(5000,1);
end

f=char(fread(fid)');fclose(fid);
% find line feed character
%i=find(f == char(10));

% find #SCAN Set(s)
ind_command = strfind(f,'#COMMAND');

fprintf('Reading data from %s... %3.0d scan set(s) found!\n',logtoopen,length(ind_command));

% initialize data cell
edd=struct('command',[],...
    'data',cell(1),...
    'motorname','',...
    'motorpos',[],...
    'motorpos_all',[],...
    'keyence',[],...
    'motor_start_end_numstep',[],...
    'motorstep',[],...
    'exp_time',[],...
    'log','');

t0 = tic;

fprintf('Start converting...');
for i = 1:length(ind_command)
    % get scan record
    if i < length(ind_command)
        snlog = f(ind_command(i):(ind_command(i+1)-1));
    else
        snlog = f(ind_command(i):(size(f,2)));
    end
    % find scans
    scanst = strfind(snlog,'#Scan');
    scanen = [scanst(2:end)-1 length(snlog)];
    % find first few line feed
    search_limit = min(1000,length(snlog));  % first 1000 characters should cover!!
    lf = find(snlog(1:search_limit) == char(10));
    % find comment
    %    comment = strfind(snlog,'#');
    command_= snlog(1:lf(1)-1);
    chk_cmd = strsplit(command_,' ');
    cmd_use = chk_cmd{2};
    
    edd(i).command = command_(11:end);
    
    variable_used = 0;  % assume no variable is used for 
    
    switch cmd_use
        case 'motor_cam'
            mname   = chk_cmd{4};
            mstart  = str2num(chk_cmd{5});
            mend    = str2num(chk_cmd{6});
            numstep = str2num(chk_cmd{7});  % number steps in command  
            exptime = str2num(chk_cmd{8});
            %%%%% if length([mstart mend mstep exptime])~=4
            % check if the scan command was finished correctly.
            nosteps = length(scanst);  % number steps recorded in file
            if nosteps ~= numstep
                warning('Scans may be interrupted!!')
            end
            % if variable is used in script, following parameters might not
            % be recorded in numerical format !!
            if isempty(mstart);  mstart = 0; variable_used = 1;end
            if isempty(mend);    mend   = 0; variable_used = 1; end
            if isempty(numstep); numstep= 0; variable_used = 1; end
            if isempty(exptime); exptime= 0; variable_used = 1; end            
            
            edd(i).log       = snlog(1:lf(3));
            edd(i).motorname = mname;
            edd(i).motor_start_end_numstep = [mstart mend numstep];
            edd(i).motorstep = (mend-mstart)/(numstep-1);
            edd(i).exp_time  = exptime;            
            edd(i).data{1}   = zeros(nosteps,8192);
            edd(i).motorpos  = zeros(1,nosteps);
            edd(i).keyence   = zeros(1,nosteps);            
            
        case 'motor_cam2'
            %%%%% To-do %%%%%
            % 1. moving motor may not be the one shown in the command.
            %%%%%%%%%%%%%%%%%
            
            mname   = chk_cmd{3};
            mstart  = str2num(chk_cmd{4});
            mend    = str2num(chk_cmd{5});
            numstep = str2num(chk_cmd{6});  % number steps in command  
            exptime = str2num(chk_cmd{7});
            %%%%% if length([mstart mend mstep exptime])~=4
            % check if the scan command was finished correctly.
            nosteps = length(scanst);  % number steps recorded in file
            if nosteps ~= numstep
                warning('Scans may be interrupted!!')
            end
            % if variable is used in script, following parameters might not
            % be recorded in numerical format !!
            if isempty(mstart);  mstart = 0; variable_used = 1;end
            if isempty(mend);    mend   = 0; variable_used = 1; end
            if isempty(numstep); numstep= 0; variable_used = 1; end
            if isempty(exptime); exptime= 0; variable_used = 1; end            
            
            edd(i).log       = snlog(1:lf(3));
            edd(i).motorname = mname;
            edd(i).motor_start_end_numstep = [mstart mend numstep];
            edd(i).motorstep = (mend-mstart)/numstep;
            edd(i).exp_time  = exptime;            
            edd(i).data      = zeros(nosteps,8192);
            edd(i).motorpos  = zeros(1,nosteps);
            edd(i).keyence   = zeros(1,nosteps);            
                        
        otherwise
            fprintf('Abort!!\n This command "%s" is not supported by this reader yet!!\n',cmd_use);
            return;
    end
    % retrieve data
    if nosteps==0
        %%% scan aborted!!
        edd(i).command = 'scan aborted!!';
    else
        for j = 1:nosteps
            tmpda = snlog(scanst(j):scanen(j));   % get data of single scan
            %comment = strfind(tmp_,'#')
            lf = find(tmpda == char(10));         % find line feed
            % header
            % header_1 = tmpda(1:lf(1)-1)
            % header_2 = tmpda(lf(1)+1:lf(2)-1)
            header_3 = tmpda(lf(2)+1:lf(3)-1);
            
            assignin('base','tmp',header_3);
            
            h3data = sscanf(header_3,'%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f %f %f %f %*f %*f %*f %*f %*f %f'); % read only keyence value
            % if keyence data exist, read it!!
            if length(h3data)==4
                edd(i).keyence(j) = h3data(4);
            end
            % save xr yr zr (keyence) position
            edd(i).motorpos_all(:,j) = h3data';
            % reassign motorposition based on moving motor
            switch edd(i).motorname
                case 'xr'
                    edd(i).motorpos(j) = h3data(1);
                case 'yr'
                    edd(i).motorpos(j) = h3data(2);
                case 'zr'
                    edd(i).motorpos(j) = h3data(3);
                case 'ksamx'   % work around for shade_oct17
                    edd(i).motorpos(j) = mstart+(j-1)*(mend-mstart)/numstep;
                otherwise
                    fprintf('\nmotor_name "%s" not supported yet!!\n',edd(i).motorname)
            end
            
            % if variable is used, calculate motor_step again from real data.
            if  variable_used
                edd(i).motor_start_end_numstep = [edd(i).motorpos(1) edd(i).motorpos(end) length(edd(i).motorpos)];
                edd(i).motorstep = round((edd(i).motorpos(end)-edd(i).motorpos(1))/(length(edd(i).motorpos)-1),3);
            end
            
            if rem(length(lf),8195) ~= 0
                warning('Scan %d may be interrupted!!',j)
            end
            da = sscanf(tmpda(lf(3)+1:end),'%f');  % read the data
            edd(i).data{1}(j,:) = da(2:2:end);
            edd(i).log = [edd(i).log tmpda(1:lf(3))];
            
            % read detector-2 data if exist
            if isdir(logtoopen)
                d2fname =fullfile(fpath,fname,sprintf('%s-%s-hv.xy',fname,padzero(j,3)));
                if exist(d2fname,'file')
                    f2id=fopen(d2fname);
                    f2=textscan(f2id,'%f %f');fclose(f2id);
                    if size(f2{2},1)< (step_counter(j)+1)*8192
                        edd(i).data{2}(j,:) = zeros(8192,1);
                    else
                        edd(i).data{2}(j,:) = f2{2}([1:8192]+step_counter(j)*8192);  % det-2 data
                    end
                    %edd(i).data3(j,:) = f2{1}([1:8192]+(i-1)*8192);  % det-1 data (save to compare with original array dump) (this seems to have higher counts than original det-1 dump)
                    %above issues has been fixed by adding a proper waiting in script (Mar.17)
                end
            end
            step_counter(j)=step_counter(j)+1;
        end
    end
end

fprintf('...done. %2.4f sec elasped.\n',toc(t0));

