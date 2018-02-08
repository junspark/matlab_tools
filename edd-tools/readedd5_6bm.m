% EDD data file reader (6BM setup)
%  read data into a single Matlab struct
%  supported configuration: after Feb.2018
%
%  Rev.1.9 (2018/02/07)
%  + read all motor information
%
%  Rev.1.8 (2018/01/18)
%  + new format to combine det-1/det-2 data in one file
%
%  Rev.1.7 (2017/10/07)
%  + read multiple motor position
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
% Copyright 2015-2018 Andrew Chuang (chuang.cp@gmail.com)
% $Revision: 1.9 $  $Date: 2018/02/07 $


function [edd]=readedd5_6bm(logtoopen)

if nargin ~= 1
    fprintf('\nUsage : [log]=readedd_6bm("FILEtoOPEN")\n');
    fprintf('\n');
    return;
end

% open and read the file into memory
if isdir(logtoopen)
    [fpath, fname, ~] = fileparts(logtoopen);
    fullname= fullfile(logtoopen,sprintf('%s.xy',fname));
    fid=fopen(fullname);if fid == -1, error('Can''t find/open the input file.'); end
    step_counter = zeros(5000,1);
else
    fid=fopen(logtoopen);if fid == -1, error('Can''t open the input file.'); end
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
    'allmotors','',...
    'motornames','',...    
    'motorname','',...
    'motorpos',[],...
    'motorpos_all',[],...
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
    lf = find(snlog(1:search_limit) == newline);  % newline = char(10)
    % find comment
    %    comment = strfind(snlog,'#');
    command_= snlog(1:lf(1)-1);
    chk_cmd = strsplit(command_,' ');
    cmd_use = chk_cmd{2};
    
    edd(i).command = command_(11:end);
    
    variable_used = 0;  % assume no variable is used to do motor_cam
    
    % name of all motors
    motor_name_all = strsplit(snlog(lf(1)+3:lf(2)-2),' ');
    % value of all motors
    motor_value_all = str2num(snlog(lf(2)+3:lf(3)-2));
    % convert all motor information into struct
    edd(i).allmotors = cell2struct(num2cell(motor_value_all),motor_name_all,2);
        
    %assignin('base','tmp',motor_all);
    
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
            if isempty(numstep)
                warning('Variable used. Can not check if scan complete or not!!')
            elseif nosteps ~= numstep
                warning('Scans may be interrupted!!')
            end
            % if variable is used in script, following parameters might not
            % be recorded in numerical format !!
            if isempty(mstart);  mstart = 0; variable_used = 1;end
            if isempty(mend);    mend   = 0; variable_used = 1; end
            if isempty(numstep); numstep= 0; variable_used = 1; end
            if isempty(exptime); exptime= 0; variable_used = 1; end            
            
            %edd(i).log       = snlog(1:lf(3));
            edd(i).log       = snlog(1:lf(1));
            edd(i).motorname = mname;
            edd(i).motor_start_end_numstep = [mstart mend numstep];
            edd(i).motorstep = (mend-mstart)/(numstep-1);
            edd(i).exp_time  = exptime;            
            edd(i).data{1}   = zeros(nosteps,8192);
            edd(i).motorpos  = zeros(1,nosteps);
%            edd(i).keyence   = zeros(1,nosteps);            
            
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
%            edd(i).keyence   = zeros(1,nosteps);            
                        
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
            lf = find(tmpda == newline);         % find line feed
            % header
            % header_1 = tmpda(1:lf(1)-1)
            % header_2 = tmpda(lf(1)+1:lf(2)-1)
            header_3 = tmpda(lf(2)+1:lf(3)-1);
            
            %assignin('base','tmp',header_3);
            
            separator = find(header_3=='=');
            % read all motor value
            h3data = str2num(header_3(separator+1:end));

            % save all motor position
            edd(i).motorpos_all(j,:) = h3data;
            
            % suppose same motor name saved for every step, save motorname only once
            if j==1
                edd(i).motornames = strsplit(header_3(2:separator-2));
                % find which motor is moving motor
                motor_no = find(~cellfun(@isempty,strfind(edd(i).motornames,edd(i).motorname)));

                if isempty(motor_no)
                    fprintf('\n"%s" position is not saved!!\n',edd(i).motorname)
                end
            end
            
            % assign motor position based on moving motor            
            if isempty(motor_no)
                edd(i).motorpos(j) = 0;
            else
                 edd(i).motorpos(j) = h3data(motor_no);
            end
            
            %motor_no = find(~cellfun(@isempty,strfind(edd(i).motornames,edd(i).motorname)))
            
            % if variable is used, calculate motor_step again from real data.
            if  variable_used
                edd(i).motor_start_end_numstep = [edd(i).motorpos(1) edd(i).motorpos(end) length(edd(i).motorpos)];
                edd(i).motorstep = round((edd(i).motorpos(end)-edd(i).motorpos(1))/(length(edd(i).motorpos)-1),3);
            end
            
            if rem(length(lf),8195) ~= 0
                warning('Scan %d may be interrupted!!',j)
            end
            da = sscanf(tmpda(lf(3)+1:end),'%f');  % read the detector data
            numcol = length(da)/8192;              % number of column
            for k = 2:numcol
                edd(i).data{k-1}(j,:) = da(k:numcol:end);     % read k-th detector data
            end
            edd(i).log = [edd(i).log tmpda(1:lf(2))];
            
            step_counter(j)=step_counter(j)+1;
        end
    end
end

fprintf('...done. %2.4f sec elasped.\n',toc(t0));

