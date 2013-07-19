function log = parseHEXRDLog(fileName)
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

fid = fopen(fileName, 'r');
if(fid == -1)
    beep;
    error('Cannot open file:\n  %s\n', fileName);
end
% Get number of grains
line = fgetl(fid);
parms = sscanf(line, 'Found %i');
nGrains = parms;
textscan(fid, '%*[^\n]', 17); % Skip 22 lines
% Loop over found grains
log(nGrains) = struct('nExpGvec',[],'nMeasGvec',[],'nMeasOnce',[], ...
    'nMeasMore',[],'meanIA',[],'U',[],'Rod',[],'Bunge',[],'Quat',[],'refl',[]);
for i = 1:nGrains
    fgetl(fid);  % Grain nnn, nPairs (Skip)
    line    = fgetl(fid);
    parms   = sscanf(line, '%d %d %d %d');
    
    log(i).nExpGvec     = parms(1);
    log(i).nMeasGvec    = parms(2);
    log(i).nMeasOnce    = parms(3);
    log(i).nMeasMore    = parms(4);
    
    line    = fgetl(fid);
    parms   = sscanf(line, '%f');  % May be more items with 0.5 and above
    
    log(i).meanIA = parms(1);
    
    C           = textscan(fid, '%f', 9);
    log(i).U    = reshape(C{1}(1:9), 3, 3)';
    
    % Skip 3 lines
    textscan(fid, '%*[^\n]', 3);
    C   = textscan(fid, '%f', 3);
    log(i).Rod  = reshape(C{1}(1:3),3,1);
    
    % Skip 1 line
    % textscan(fid, '%*[^\n]', 1)
    C               = textscan(fid, '%f', 3);
    log(i).Bunge    = reshape(C{1}(1:3),3,1);
    
    C           = textscan(fid, '%f', 4);
    log(i).Quat = reshape(C{1}(1:4),4,1);
    
    % gvec is the whole 22 x nMeasGvec array of columns
    C           = textscan(fid, '%f', 22*log(i).nMeasGvec);
    log(i).refl = reshape(C{1},22,log(i).nMeasGvec)';
    textscan(fid, '%*[^\n]', 1); % Skip 1 lines
end
fclose(fid);
