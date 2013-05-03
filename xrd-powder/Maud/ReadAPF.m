function apf_data = ReadAPF(fname)
% ReadAPF - reads APF file from MAUD
%
%   USAGE:
%
%   ReadAPF(fname)
%
%   INPUT:
%
%   fname
%           name of the APF file
%
%   OUTPUT:
%
%   apf_data
%			is the pole figure data extracted from 
%			x-ray diffraction spectrum after MAUD analysis	
%

fid     = fopen(fname, 'r');

linedata	= fgetl(fid);
idx         = strfind(linedata,' ');
idx         = idx(5) + 1;
apf_data.phase_name = linedata(idx:end);

linedata	= fgetl(fid);
idx         = strfind(linedata,' ');
apf_data.num_hkl    = str2double(linedata(1:idx(1)));

for i = 1:1:apf_data.num_hkl
    linedata	= fgetl(fid);
    idx = strfind(linedata,'   ');
    h   = str2double(linedata(1:idx(1)));
    k   = str2double(linedata(idx(1):idx(2)));
    l   = str2double(linedata(idx(2):idx(3)));
    apf_data.pf(i).hkl  = [h k l];
    
    linedata    = fgetl(fid);
    idx = strfind(linedata, ' <-  MEPSUM');
    apf_data.pf(i).num_pts  = str2double(linedata(1:idx(1)));
    pfdata  = zeros(apf_data.pf(i).num_pts, 11);
    for j = 1:1:apf_data.pf(i).num_pts
        linedata    = fgetl(fid);
        pfdata(j,:) = str2num(linedata);
    end
    apf_data.pf(i).polar_angle  = pfdata(:,1);
    apf_data.pf(i).azi_angle    = pfdata(:,2);
    apf_data.pf(i).value_from_exp   = pfdata(:,3);
    apf_data.pf(i).value_from_odf   = pfdata(:,4);
    apf_data.pf(i).prog_number      = pfdata(:,5);
    apf_data.pf(i).weight           = pfdata(:,6);
    
    apf_data.pf(i).omega        = pfdata(:,7);
    apf_data.pf(i).chi          = pfdata(:,8);
    apf_data.pf(i).phi          = pfdata(:,9);
    apf_data.pf(i).eta          = pfdata(:,10);
    apf_data.pf(i).bank_num     = pfdata(:,11);
end
fclose(fid);
