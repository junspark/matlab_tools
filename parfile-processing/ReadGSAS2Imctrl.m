function imctrl = ReadGSAS2Imctrl(fimctrl)

% clear all
% close all
% clc
% 
% fclose all;
% 
% fimctrl = 'mli_dec19_tif_am316_ss_id105_1.imctrl';

fid = fopen(fimctrl, 'r');
while ~feof(fid)
    lindata = fgetl(fid);
    
    if contains(lindata, 'LRazimuth')
        idx1    = strfind(lindata, '[');
        idx2    = strfind(lindata, ']');
        imctrl.azm_range    = sscanf(lindata(idx1+1:idx2-1), '%f, %f');
    elseif contains(lindata, 'azmthOff')
        idx1    = strfind(lindata, ':');
        imctrl.azm_offset   = sscanf(lindata(idx1+1:end), '%f');
    elseif contains(lindata, 'outAzimuths')
        idx1    = strfind(lindata, ':');
        imctrl.azm_steps    = sscanf(lindata(idx1+1:end), '%d');
    elseif contains(lindata, 'IOtth')
        idx1    = strfind(lindata, '[');
        idx2    = strfind(lindata, ']');
        imctrl.tth_range    = sscanf(lindata(idx1+1:idx2-1), '%f, %f');
    elseif contains(lindata, 'outChannels')
        idx1    = strfind(lindata, ':');
        imctrl.tth_steps    = sscanf(lindata(idx1+1:end), '%d');
    end
end
fclose(fid);
