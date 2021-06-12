function pfname_polimg = GSAS2_MakePolImgFromChi(pfname_imctrl, pname_chi, froot_chi)
% GSAS2_MakePolImgFromChi - Makes PolImg file from GSAS2 PWDR CHI files
%
%   INPUT:
%
%   pfname_imctrl
%       full pathname to the imctrl file used for integration.
%   
%   pname_chi
%       pathname to where the chi files exist.
%
%   froot_chi
%       file name root of rht chi files.
%
%   OUTPUT:
%
%   pfname_polimg
%       full pathname to where the synthetic polimg file was saved.

% clear all
% close all
% clc

%%% INPUTS
% pfname_imctrl   = 'D:\w\mpe_apr21_analysis\park_sae_phi0_2\mpe_apr21_park_sae_ge3_80keV_2278mm_phi0.imctrl';
% pname_chi       = 'D:\w\mpe_apr21_analysis\park_sae_phi0_2\chi';
% froot_chi       = 'park_sae_phi0_2_001749.ge3.sum';
%%%

%%% GENERATE METADATA
imctrl  = GSAS2_Read_imctrl(pfname_imctrl);

% dazm        = imctrl.dazm;
azm_grid	= imctrl.derived_azm_grid;

%%%% cakeprms
cakeprms.fastint	= 2;
cakeprms.bins       = [imctrl.outAzimuths imctrl.outChannels nan];
cakeprms.origin     = imctrl.center;
cakeprms.sector     = [azm_grid(1) azm_grid(end) imctrl.IOtth(1) imctrl.IOtth(2)];
cakeprms.azim       = azm_grid;

%%%% instr
instr.wavelength    = imctrl.wavelength;
instr.energy        = Angstrom2keV(instr.wavelength);
instr.distance      = imctrl.distance;
instr.centers       = [nan nan];
instr.gammaX        = nan;
instr.gammaY        = nan;
instr.detsizeHorz	= nan;
instr.detsizeVert   = nan;
instr.pixelsizeHorz = nan;
instr.pixelsizeVert	= nan;
instr.numpixelsHorz	= nan;
instr.numpixelsVert	= nan;
instr.imrotation    = nan;
instr.dettype       = nan;
instr.detpars       = nan;
instr.omega         = nan;
instr.chi           = nan;
instr.prrot         = nan;

for iii = 1:1:cakeprms.bins(1)
    % fname_chi   = sprintf('%s_Azm=_%d.chi', froot_chi, fix(cakeprms.azim(iii)));
    fname_chi   = sprintf('%s_A%d.chi', froot_chi, fix(cakeprms.azim(iii)));
    pfname_chi  = fullfile(pname_chi, fname_chi);
    
    data_chi    = GSAS2_Read_PWDR_CHI(pfname_chi);
    
    mapped_tth_for_intensity(iii,:) = data_chi(:,1);
    intensity_in_tth_grid(iii,:)    = data_chi(:,2);
end

polimg.azimuth                  = cakeprms.azim;
polimg.radius                   = instr.distance.*tand(mapped_tth_for_intensity);
polimg.intensity                = intensity_in_tth_grid;
polimg.mapped_tth_for_intensity = mapped_tth_for_intensity;
polimg.tth_grid                 = mean(mapped_tth_for_intensity,1);
polimg.d_grid                   = (instr.wavelength/2)./sind(polimg.tth_grid./2);
polimg.intensity_in_tth_grid    = intensity_in_tth_grid;

fname_polimg    = sprintf('%s.polimg.mat', froot_chi);
pfname_polimg   = fullfile(pname_chi, fname_polimg);
save(pfname_polimg, 'polimg', 'instr', 'cakeprms');