function log = parseHEXRDLog(pfname)
% parseHEXRDLog Parse a Grainspotter log-like HEXRD log file.
%
%   log = parseGrainSpotterLog(fileName) reads the Grainspotter log file
%   with the name fileName and returns the information in an array of
%   structures with fields:
%       nExpGvec = Number of expected G vectors
%       nMeasGvec = Number of measured G vectors
%       nMeasOnce = Number of G vectors measured once
%       nMeasMore = Number of G vectors measured more than once
%       meanIA = Average internal angle between prediced and measured
%       U = 3x3 Orientation matrix
%       gvec = G vector table
%       hkl = 3 hkl values
%
%   The columns of the log file are:
%     n gvector_id peak_id  h k l  h_pred k_pred l_pred  dh dk dl
%     tth_meas tth_pred dtth  omega_meas omega_pred domega
%     eta_meas  eta_pred deta  IA
%
%   Example:
%     log = parseGrainSpotterLog('simul.log');

[fid, msg]  = fopen(pfname, 'r');
if(fid == -1)
    beep;
    error('Cannot open file:\n  %s\n', fileName);
end

nGrains = 0;
while ~feof(fid)
    line = fgetl(fid);
    if strfind(line, '#### grain') == 1
        nGrains  = nGrains + 1;
    end
end

% Rewind file
fseek(fid, 0, -1);

% Loop over found grains
log(nGrains) = struct(...
    'Quat',[],'R',[],'V',[],'Esam',[],'Ecry',[],'F',[], ...
    'lattprms',[], 'COM', [], 'ReflectionTable', [], 'Completeness', []);

while ~feof(fid)
    line = fgetl(fid);
    if strfind(line, '#### grain') == 1
        grnum   = str2double(line(11:end)) + 1;
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 1;
        log(grnum).Quat = str2num(line(idx:end-1));
        
        %%% LINE SKIP
        fgetl(fid);
        
        R       = zeros(3,3);
        line    = fgetl(fid);
        idx     = strfind(line, '[[') + 2;
        R(1,:)  = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx 	= strfind(line, '[') + 1;
        R(2,:)  = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 2;
        R(3,:)  = str2num(line(idx:end-2));
        log(grnum).R    = R;
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        
        V       = zeros(3,3);
        line    = fgetl(fid);
        idx     = strfind(line, '[[') + 2;
        V(1,:)  = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 1;
        V(2,:)  = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 2;
        V(3,:)  = str2num(line(idx:end-2));
        log(grnum).V    = V;
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        
        Esam    = zeros(3,3);
        line    = fgetl(fid);
        idx     = strfind(line, '[[') + 2;
        Esam(1,:)   = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 1;
        Esam(2,:)   = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 2;
        Esam(3,:)   = str2num(line(idx:end-2)); 
        log(grnum).Esam = Esam;
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        
        Ecry    = zeros(3,3);
        line    = fgetl(fid);
        idx     = strfind(line, '[[') + 2;
        Ecry(1,:)   = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 1;
        Ecry(2,:)   = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 2;
        Ecry(3,:)   = str2num(line(idx:end-2)); 
        log(grnum).Ecry = Ecry;
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        
        F       = zeros(3,3);
        line    = fgetl(fid);
        idx     = strfind(line, '[[') + 2;
        F(1,:)  = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 1;
        F(2,:)  = str2num(line(idx:end-2));
        
        line    = fgetl(fid);
        idx     = strfind(line, '[') + 2;
        F(3,:)  = str2num(line(idx:end-2)); 
        log(grnum).F    = F;
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        
        line    = fgetl(fid);
        log(grnum).lattprms = str2num(line(2:end));
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        
        line    = fgetl(fid);
        idx     = strfind(line, '(') + 1;
        log(grnum).COM  = str2num(line(idx:end-1));
        
        %%% LINE SKIP
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        
        EndOfTable  = 0;
        ReflectionTable = [];
        while ~EndOfTable
            line    = fgetl(fid);
            ReflectionTable = [ReflectionTable; str2num(line)];
            if strfind(line, '#')
                EndOfTable  = 1;
            end
        end
        log(grnum).ReflectionTable  = ReflectionTable;
        
        %%% LINE SKIP
        line    = fgetl(fid);
        idx     = strfind(line, ':');
        log(grnum).Completeness = str2double(line(idx+1:end-1));
    end
end
fclose(fid);