clear all
close all
clc

addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%% DATA REDUCTION FLAGS
Analysis_Options.PkFitOptimizationOptions   = optimset(...
    'MaxIter', 5e5, ...
    'MaxFunEvals',3e5);

Analysis_Options.make_polimg    = 0;
Analysis_Options.save_polimg    = 0;
Analysis_Options.fits_spectra   = 1;
Analysis_Options.save_fits      = 1;
Analysis_Options.find_instrpars = 0;
Analysis_Options.save_instrpars = 0;
Analysis_Options.find_detpars	= 0;
Analysis_Options.generateESG    = 0;

%%% INPUT PARAMETERS
XRDIMAGE.ExpID              = 'mli_dec19';
XRDIMAGE.MetaDataFile       = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\mli_nov19_Tomo.par';
XRDIMAGE.Image.pname_chi    = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\data\dexela\corrected';
XRDIMAGE.Image.pname        = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\data\dexela\corrected';
XRDIMAGE.Image.fbase        = 'am316_ss_id105_1';
XRDIMAGE.Image.fnumber      = 324:326;
XRDIMAGE.Image.ref_fnumber  = 324;   %%% REFERENCE INDEX
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.corrected    = 1;   % 0 - UNCORRECTED, 1 = SUM, 2 = COR32
XRDIMAGE.Image.dettype      = 6;   % 3 = GE3, 5 = GE5, 1234 = HYDRA1234, 6 = DEX IN C

%%% FITTING PARAMETER
XRDIMAGE.pkListFile         = 'C:\Users\parkjs\Documents\MATLAB\work\mli_nov19_analysis\waxs\am316_ss_id105_1.pkslst';

%%% REPROCESSED FILE NAME
fname_reprocessed   = sprintf('%s.reprocessed.mat', XRDIMAGE.Image.fbase);

%%% ANALYSIS STARTS HERE
metadata    = ReadSpecParFile(XRDIMAGE.MetaDataFile, 'Version', 'mli_nov19_c');

%%% READ PEAK LIST FROM GSAS2
fid     = fopen(XRDIMAGE.pkListFile, 'r');
pklist  = [];
while ~feof(fid)
    linedata    = fgetl(fid);
    switch linedata(1)
        case '['
            pklist  = [pklist; eval(linedata)];
    end
end
fclose(fid);
numpks  = size(pklist ,1);

switch XRDIMAGE.Image.dettype
    case 6
        disp(sprintf('%d = DEXELA', XRDIMAGE.Image.dettype));
        %%% FIRST PASS - SIMPLY LOAD / REARRANGE DATA
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            
            %%% eta x 1
            eta{iii}    = [];
            ef{iii}     = [];
            amp{iii}    = [];
            fwhm{iii}   = [];
            mix{iii}    = [];
            rho{iii}    = [];
            rsn{iii}    = [];
            rwp{iii}    = [];
            
            if XRDIMAGE.Image.fnumber(iii) == XRDIMAGE.Image.ref_fnumber
                amp0    = [];
                fwhm0   = [];
                mix0    = [];
                rho0    = [];
            end
            
            fext   = 'tif';
            
            fimctrl = sprintf('%s_%s_%s.imctrl', XRDIMAGE.ExpID, fext, XRDIMAGE.Image.fbase);
            imctrl  = ReadGSAS2Imctrl(fimctrl);
            
            dazm        = (imctrl.azm_range(2) - imctrl.azm_range(1))/imctrl.azm_steps;
            azm_grid    = linspace(imctrl.azm_range(1)+dazm/2, imctrl.azm_range(2)-dazm/2, imctrl.azm_steps);
            
            idx_azm             = azm_grid >= 360;
            azm_grid(idx_azm)   = azm_grid(idx_azm) - 360;
            switch XRDIMAGE.Image.corrected
                case 1
                    fname_pattern   = sprintf('%s_bkg_corrected_%s', froot, fext);
                    fname_pkfit     = sprintf('%s.pkfit.mat',fname_pattern);
                    pfname_pkfit    = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fname_pkfit);
                    
                    fprintf('processing %s ...\n', fname_pkfit)
                    switch isfile(pfname_pkfit)
                        case true
                            %%% START LOADING PKFITS
                            pkfit_data = load(pfname_pkfit);

                            if XRDIMAGE.Image.fnumber(iii) == XRDIMAGE.Image.ref_fnumber
                                amp0    = [amp0; pkfit_data.pkfit_amp];
                                fwhm0   = [fwhm0; pkfit_data.pkfit_fwhm];
                                mix0    = [mix0; pkfit_data.pkfit_mix];
                                rho0    = [rho0; pkfit_data.pkfit_rho];
                            end

                            eta{iii}    = [eta{iii}; pkfit_data.pkfit_eta];

                            ef{iii}     = [ef{iii}; pkfit_data.pkfit_ef];
                            amp{iii}    = [amp{iii}; pkfit_data.pkfit_amp];
                            fwhm{iii}   = [fwhm{iii}; pkfit_data.pkfit_fwhm];
                            mix{iii}    = [mix{iii}; pkfit_data.pkfit_mix];
                            rho{iii}    = [rho{iii}; pkfit_data.pkfit_rho];
                            rsn{iii}    = [rsn{iii}; pkfit_data.pkfit_rsn];
                            rwp{iii}    = [rwp{iii}; pkfit_data.pkfit_rwp];

                            rX(iii) = pkfit_data.pkfit_rX;
                            rY(iii) = pkfit_data.pkfit_rY;
                            rZ(iii) = pkfit_data.pkfit_rZ;

%                             samP(iii)   = pkfit_data.pkfit_samP;
%                             samX(iii)   = pkfit_data.pkfit_samX;
%                             samX2(iii)  = pkfit_data.pkfit_samX2;
%                             samY(iii)   = pkfit_data.pkfit_samY;
%                             samY2(iii)  = pkfit_data.pkfit_samY2;
%                             samZ(iii)   = pkfit_data.pkfit_samZ;
%                             samZ2(iii)  = pkfit_data.pkfit_samZ2;
% 
%                             enc(1,iii)  = pkfit_data.pkfit_enc1;
%                             enc(2,iii)  = pkfit_data.pkfit_enc2;
%                             enc(3,iii)  = pkfit_data.pkfit_enc3;
%                             enc(4,iii)  = pkfit_data.pkfit_enc4;
%                             enc(5,iii)  = pkfit_data.pkfit_enc5;ef{iii, jjj}    = pkfit_data.pkfit_ef;
%                             enc(6,iii)  = pkfit_data.pkfit_enc6;
%                             enc(7,iii)  = pkfit_data.pkfit_enc7;
%                             enc(8,iii)  = pkfit_data.pkfit_enc8;
%                             enc(9,iii)  = pkfit_data.pkfit_enc9;
%                             enc(10,iii) = pkfit_data.pkfit_enc10;
% 
%                             ev(1,iii)   = pkfit_data.pkfit_ev1;
%                             ev(2,iii)   = pkfit_data.pkfit_ev2;
%                             ev(3,iii)   = pkfit_data.pkfit_ev3;
%                             ev(4,iii)   = pkfit_data.pkfit_ev4;
%                             ev(5,iii)   = pkfit_data.pkfit_ev5;
%                             ev(6,iii)   = pkfit_data.pkfit_ev6;
%                             ev(7,iii)   = pkfit_data.pkfit_ev7;
%                             ev(8,iii)   = pkfit_data.pkfit_ev8;
%                             ev(9,iii)   = pkfit_data.pkfit_ev9;
%                             ev(10,iii)  = pkfit_data.pkfit_ev10;
                        otherwise
                            eta{iii}    = [eta{iii}; nan(azm_steps, 1)];

                            ef{iii}     = [ef{iii}; nan(azm_steps, numpks)];
                            amp{iii}    = [amp{iii}; nan(azm_steps, numpks)];
                            fwhm{iii}   = [fwhm{iii}; nan(azm_steps, numpks)];
                            mix{iii}    = [mix{iii}; nan(azm_steps, numpks)];
                            rho{iii}    = [rho{iii}; nan(azm_steps, numpks)];
                            rsn{iii}    = [rsn{iii}; nan(azm_steps, numpks)];
                            rwp{iii}    = [rwp{iii}; nan(azm_steps, numpks)];

                            rX(iii) = nan;
                            rY(iii) = nan;
                            rZ(iii) = nan;

%                             samP(iii)   = nan;
%                             samX(iii)   = nan;
%                             samX2(iii)  = nan;
%                             samY(iii)   = nan;
%                             samY2(iii)  = nan;
%                             samZ(iii)   = nan;
%                             samZ2(iii)  = nan;
% 
%                             enc(1,iii)  = nan;
%                             enc(2,iii)  = nan;
%                             enc(3,iii)  = nan;
%                             enc(4,iii)  = nan;
%                             enc(5,iii)  = nan;
%                             enc(6,iii)  = nan;
%                             enc(7,iii)  = nan;
%                             enc(8,iii)  = nan;
%                             enc(9,iii)  = nan;
%                             enc(10,iii) = nan;
% 
%                             ev(1,iii)   = nan;
%                             ev(2,iii)   = nan;
%                             ev(3,iii)   = nan;
%                             ev(4,iii)   = nan;
%                             ev(5,iii)   = nan;
%                             ev(6,iii)   = nan;
%                             ev(7,iii)   = nan;
%                             ev(8,iii)   = nan;
%                             ev(9,iii)   = nan;
%                             ev(10,iii)  = nan;
                    end
                otherwise
                    warning('XRDIMAGE.Image.corrected = %d is not a valid option', XRDIMAGE.Image.corrected);
                    return
            end
        end
        
        %%% SECOND PASS
        %%% COMPUTE STRAINS
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            fprintf('computing lattice strain for %s .... \n', froot)
            lattice_strain{iii} = sind(rho0./2)./sind(rho{iii}./2) - 1;
            
            %%% SCATTERING VECTOR CALCULATION HERE
            % q{iii}              =
        end
    case 1234
        sprintf('%d = HYDRA1234', XRDIMAGE.Image.dettype);
        
        %%% FIRST PASS - SIMPLY LOAD / REARRANGE DATA
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            
            %%% eta x 1
            eta{iii}    = [];
            ef{iii}     = [];
            amp{iii}    = [];
            fwhm{iii}   = [];
            mix{iii}    = [];
            rho{iii}    = [];
            rsn{iii}    = [];
            rwp{iii}    = [];
            
            if XRDIMAGE.Image.fnumber(iii) == XRDIMAGE.Image.ref_fnumber
                amp0    = [];
                fwhm0   = [];
                mix0    = [];
                rho0    = [];
            end
            
            for jjj = 1:1:4
                fext    = sprintf('ge%d', jjj);
                fimctrl = sprintf('%s_%s_%s.imctrl', XRDIMAGE.ExpID, fext, XRDIMAGE.Image.fbase);
                
                imctrl  = readtable(fimctrl, 'FileType', 'text', 'ReadVariableNames', false, ...
                    'delimiter', ':');
                azm_steps   = str2double(cell2mat(imctrl.Var2(find(ismember(imctrl.Var1, 'outAzimuths')))));
                switch XRDIMAGE.Image.corrected
                    case 1
                        fname_pattern   = sprintf('%s.%s.sum', froot, fext);
                        
                        fname_pkfit    = sprintf('%s.pkfit.mat',fname_pattern);
                        pfname_pkfit    = fullfile(XRDIMAGE.Image.pname_chi, XRDIMAGE.Image.fbase, fname_pkfit);
                        
                        fprintf('processing %s ...\n', fname_pkfit)
                        switch isfile(pfname_pkfit)
                            case true
                                %%% START LOADING PKFITS
                                pkfit_data = load(pfname_pkfit);
                                % pkfit_data
                                
                                if XRDIMAGE.Image.fnumber(iii) == XRDIMAGE.Image.ref_fnumber
                                    amp0    = [amp0; pkfit_data.pkfit_amp];
                                    fwhm0   = [fwhm0; pkfit_data.pkfit_fwhm];
                                    mix0    = [mix0; pkfit_data.pkfit_mix];
                                    rho0    = [rho0; pkfit_data.pkfit_rho];
                                end
                                
                                eta{iii}    = [eta{iii}; pkfit_data.pkfit_eta];
                                
                                ef{iii}     = [ef{iii}; pkfit_data.pkfit_ef];
                                amp{iii}    = [amp{iii}; pkfit_data.pkfit_amp];
                                fwhm{iii}   = [fwhm{iii}; pkfit_data.pkfit_fwhm];
                                mix{iii}    = [mix{iii}; pkfit_data.pkfit_mix];
                                rho{iii}    = [rho{iii}; pkfit_data.pkfit_rho];
                                rsn{iii}    = [rsn{iii}; pkfit_data.pkfit_rsn];
                                rwp{iii}    = [rwp{iii}; pkfit_data.pkfit_rwp];
                                
                                rX(iii) = pkfit_data.pkfit_rX;
                                rY(iii) = pkfit_data.pkfit_rY;
                                rZ(iii) = pkfit_data.pkfit_rZ;
                                
                                samP(iii)   = pkfit_data.pkfit_samP;
                                samX(iii)   = pkfit_data.pkfit_samX;
                                samX2(iii)  = pkfit_data.pkfit_samX2;
                                samY(iii)   = pkfit_data.pkfit_samY;
                                samY2(iii)  = pkfit_data.pkfit_samY2;
                                samZ(iii)   = pkfit_data.pkfit_samZ;
                                samZ2(iii)  = pkfit_data.pkfit_samZ2;
                                
                                enc(1,iii)  = pkfit_data.pkfit_enc1;
                                enc(2,iii)  = pkfit_data.pkfit_enc2;
                                enc(3,iii)  = pkfit_data.pkfit_enc3;
                                enc(4,iii)  = pkfit_data.pkfit_enc4;
                                enc(5,iii)  = pkfit_data.pkfit_enc5;ef{iii, jjj}    = pkfit_data.pkfit_ef;
                                enc(6,iii)  = pkfit_data.pkfit_enc6;
                                enc(7,iii)  = pkfit_data.pkfit_enc7;
                                enc(8,iii)  = pkfit_data.pkfit_enc8;
                                enc(9,iii)  = pkfit_data.pkfit_enc9;
                                enc(10,iii) = pkfit_data.pkfit_enc10;
                                
                                ev(1,iii)   = pkfit_data.pkfit_ev1;
                                ev(2,iii)   = pkfit_data.pkfit_ev2;
                                ev(3,iii)   = pkfit_data.pkfit_ev3;
                                ev(4,iii)   = pkfit_data.pkfit_ev4;
                                ev(5,iii)   = pkfit_data.pkfit_ev5;
                                ev(6,iii)   = pkfit_data.pkfit_ev6;
                                ev(7,iii)   = pkfit_data.pkfit_ev7;
                                ev(8,iii)   = pkfit_data.pkfit_ev8;
                                ev(9,iii)   = pkfit_data.pkfit_ev9;
                                ev(10,iii)  = pkfit_data.pkfit_ev10;
                            otherwise
                                eta{iii}    = [eta{iii}; nan(azm_steps, 1)];
                                
                                ef{iii}     = [ef{iii}; nan(azm_steps, numpks)];
                                amp{iii}    = [amp{iii}; nan(azm_steps, numpks)];
                                fwhm{iii}   = [fwhm{iii}; nan(azm_steps, numpks)];
                                mix{iii}    = [mix{iii}; nan(azm_steps, numpks)];
                                rho{iii}    = [rho{iii}; nan(azm_steps, numpks)];
                                rsn{iii}    = [rsn{iii}; nan(azm_steps, numpks)];
                                rwp{iii}    = [rwp{iii}; nan(azm_steps, numpks)];
                                
                                rX(iii) = nan;
                                rY(iii) = nan;
                                rZ(iii) = nan;
                                
                                samP(iii)   = nan;
                                samX(iii)   = nan;
                                samX2(iii)  = nan;
                                samY(iii)   = nan;
                                samY2(iii)  = nan;
                                samZ(iii)   = nan;
                                samZ2(iii)  = nan;
                                
                                enc(1,iii)  = nan;
                                enc(2,iii)  = nan;
                                enc(3,iii)  = nan;
                                enc(4,iii)  = nan;
                                enc(5,iii)  = nan;
                                enc(6,iii)  = nan;
                                enc(7,iii)  = nan;
                                enc(8,iii)  = nan;
                                enc(9,iii)  = nan;
                                enc(10,iii) = nan;
                                
                                ev(1,iii)   = nan;
                                ev(2,iii)   = nan;
                                ev(3,iii)   = nan;
                                ev(4,iii)   = nan;
                                ev(5,iii)   = nan;
                                ev(6,iii)   = nan;
                                ev(7,iii)   = nan;
                                ev(8,iii)   = nan;
                                ev(9,iii)   = nan;
                                ev(10,iii)  = nan;
                        end
                    otherwise
                        warning('XRDIMAGE.Image.corrected = %d is not a valid option', XRDIMAGE.Image.corrected);
                        return
                end
            end
        end
        
        %%% SECOND PASS
        %%% MAKE HYDRA INTO ONE
        %%% COMPUTE STRAINS
        for iii = 1:1:length(XRDIMAGE.Image.fnumber)
            froot   = sprintf('%s_%06d', XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(iii));
            fprintf('computing lattice strain for %s .... \n', froot)
            lattice_strain{iii} = sind(rho0./2)./sind(rho{iii}./2) - 1;
            
            %%% SCATTERING VECTOR CALCULATION HERE
            % q{iii}              =
        end
    otherwise
        warning('%d is not a valid option', XRDIMAGE.Image.dettype);
end

switch isfile(fname_reprocessed)
    case true
        disp('###########################')
        fprintf('Reprocessed file already exists %s\n', fname_reprocessed)
        fprintf('Rename the old file if it needs to be preserved.\n')
        fprintf('When ready, press any key to continue.\n')
        pause
    otherwise
        disp('###########################')
        fprintf('Saving reprocessed file %s\n', fname_reprocessed)
        save(fname_reprocessed)
end
