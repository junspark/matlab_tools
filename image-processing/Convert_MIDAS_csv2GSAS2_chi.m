function [] = Convert_MIDAS_csv2GSAS2_chi(pfname_MIDAS_csv, pfname_GSAS2_chi, fname_GSAS2_chi, each_azim_out)
% Convert_MIDAS_csv2GSAS2_chi - Converts MIDAS CSV TO GSAS2 CHI
%
%   USAGE:
%
%   Convert_MIDAS_csv2GSAS2_chi(pfname_MIDAS_csv, pfname_GSAS2_chi, fname_GSAS2_chi, each_azim_out)
%
%   INPUT:
%
%   pfname_MIDAS_csv        full path where MIDAS CSV file is located.
%
%   pfname_GSAS2_chi        full path where output chi file will be saved.
%
%   fname_GSAS2_chi         header information (filename of the chi file).
%
%   each_azim_out           output each azimuth? 

switch each_azim_out 
    case 1
        disp('not implemented yet')
    case 0
        data_MIDAS_csv  = load(pfname_MIDAS_csv);
        radius          = data_MIDAS_csv(:,1);
        tth             = data_MIDAS_csv(:,2);
        eta             = data_MIDAS_csv(:,3);
        intensity       = data_MIDAS_csv(:,4);
        
        nazim   = length(unique(eta));
        ntth    = length(unique(radius));
        
        radius      = mean(reshape(radius, nazim, ntth)', 2);
        tth         = mean(reshape(tth, nazim, ntth)', 2);
        eta         = mean(reshape(eta, nazim, ntth)', 2);
        intensity   = sum(reshape(intensity, nazim, ntth)', 2);
        
        % plot(tth, intensity, '.-')
        
        fid = fopen(pfname_GSAS2_chi, 'w');
        fprintf(fid, '%s Azm= 0.0\n', fname_GSAS2_chi);
        fprintf(fid, '2-Theta Angle (Degrees)\n');
        fprintf(fid, 'Intensity\n');
        fprintf(fid, '       %d\n', ntth);
        fprintf(fid, ' %f %f\n', [tth intensity]');
        fclose(fid);
        
    otherwise
        disp('choice is true or false (1 / 0)')
end
            

