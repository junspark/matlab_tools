function BuildESG(pfname, polImg, distance, wavelength, omega, chi, cakeParms)
% BuildESG - builds ESG file for MAUD input
%
%   USAGE:
%
%   BuildESG(pfname, polImg, distance, wavelength, omega, chi, cakeParms)
%
%   INPUT:
%
%   pfname
%           full file name of the output esg file
% 
%   polImg
%           is the caked data generated using instr and cakeParms
% 
%   distance
%           is the sample to detector distance
%
%   wavelength
%           is the sample to detector distance
%
%   omega
%           is the rotation about the Y axis in APS
%
%   chi
%           is the rotation about the Z axis in APS (rarely used)
%
%   cakeParms
%           is the structure array that defines the caking parameters
%
%   OUTPUT:  
%           
%           none
%

instr.distance      = distance;
instr.wavelength    = wavelength;
instr.omega         = omega;
instr.chi           = chi;

r   = tand(polImg.tth_grid).*instr.distance;
I   = polImg.intensity_in_tth_grid;

disp(['writing ', pfname])
fid = fopen(pfname,'w');
A   = {'_pd_block_id noTitle|#0';
    '';
    '_diffrn_detector Image Plate';
    '_diffrn_detector_type Image Plate';
    '_pd_meas_step_count_time ?';
    '_diffrn_measurement_method ?';
    '_diffrn_measurement_distance_unit mm';
    ['_pd_instr_dist_spec/detc ' num2str(instr.distance)];
    ['_diffrn_radiation_wavelength ' num2str(instr.wavelength)];
    '_diffrn_source_target ?';
    '_diffrn_source_power ?';
    '_diffrn_source_current ?';
    ['_pd_meas_angle_omega ' num2str(instr.omega)];
    ['_pd_meas_angle_chi ' num2str(instr.chi)];
    '_pd_meas_angle_phi 0.0';
    '_riet_par_spec_displac_x 0';
    '_riet_par_spec_displac_y 0';
    '_riet_par_spec_displac_z 0';
    '_riet_meas_datafile_calibrated false';
    '_pd_meas_angle_eta 0.0';
    '';
    'loop_';
    '_pd_proc_2theta_corrected';
    '_pd_calc_intensity_total';};

for yyy=1:size(A,1);
    fprintf(fid, '%s', A{yyy});
    fprintf(fid, '\n');
end

for xxx = 1:1:cakeParms.bins(1)
    for yyy = 1:1:size(r,2)
        fprintf(fid, '%d %d', [r(1,size(r,2)-yyy+1) I(xxx,size(r,2)-yyy+1)]);
        fprintf(fid, '\n');
    end
    
    if xxx < cakeParms.bins(1)
        A={'';
            ['_pd_block_id noTitle|#' num2str(xxx)];
            '';
            ['_pd_meas_angle_eta ' num2str(cakeParms.azim(xxx+1))];
            '';
            'loop_';
            '_pd_proc_2theta_corrected';
            '_pd_calc_intensity_total';};
        for yyy = 1:1:size(A,1);
            fprintf(fid, '%s', A{yyy});
            fprintf(fid, '\n');
        end
    end
end
fclose(fid);
