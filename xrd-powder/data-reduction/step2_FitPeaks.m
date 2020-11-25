clear all
close all
clc


addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%% CALIBRATION FILES
XRDIMAGE.Calib.pname        = '/home/beams/S1IDUSER/mnt/orthros2/mpe1_oct20_bc';
XRDIMAGE.Calib.froot        = 'LaB6_1s';
XRDIMAGE.Calib.fnumber      = 6:7; % 4116 / 4117
XRDIMAGE.Calib.numdigs      = 6;
XRDIMAGE.Calib.fext         = 'ge3';
XRDIMAGE.Calib.corrected    = 'sum';

%%% SAMPLE FILES
XRDIMAGE.Image.pname        = '/home/beams/S1IDUSER/mnt/orthros2/mpe1_oct20_bc';
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.fext         = 'ge3';        %%% CHECK METADATA SEARCH
XRDIMAGE.Image.corrected    = 'cor32';
XRDIMAGE.Image.IsHydra      = 0;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% METADATA FILES
XRDIMAGE.Image.per_scan_metadata    = '/home/beams/S1IDUSER/new_data/mpe1_oct20/saxs_waxs_fmt_fastpar.par';
XRDIMAGE.Image.per_frame_metadata   = '/home/beams/S1IDUSER/new_data/mpe1_oct20/mpe1_oct20_FF.par';

numimg  = length(XRDIMAGE.Calib.fnumber);
for iii = 1:1:numimg
    fname   = sprintf('%s_%06d.%s.%s', ...
        XRDIMAGE.Calib.froot, XRDIMAGE.Calib.fnumber(iii), XRDIMAGE.Calib.fext, XRDIMAGE.Calib.corrected);
    pfname{iii, 1}  = fullfile(XRDIMAGE.Calib.pname, ...
        XRDIMAGE.Calib.froot, XRDIMAGE.Calib.fext, fname);
end

%%% INSTRUMENT PARAMETERS GETS LOADED TO DO SOME INITIAL CALCULATION
Instr.energy        = 0;
Instr.wavelength	= 0;
Instr.distance      = 0;
Instr.centers       = [0, 0];
Instr.gammaX        = 0;
Instr.gammaY        = 0;
Instr.detsizeHorz	= 0;
Instr.detsizeVert   = 0;
Instr.pixelsizeHorz	= 0;
Instr.pixelsizeVert	= 0;
Instr.numpixelsHorz	= 0;
Instr.numpixelsVert	= 0;
Instr.imrotation    = 0;
Instr.dettype       = '2a';
Instr.detpars       = [0 0 0 0 0 0];
for iii = 1:1:length(XRDIMAGE.Calib.fnumber)
    pfname_instr    = [pfname{iii,1}, '.instr.mat'];
    Instri  = load(pfname_instr);
    Instri  = Instri.Instr;
    
    Instr.energy        = Instr.energy + Instri.energy;
    Instr.wavelength    = Instr.wavelength + Instri.wavelength;
    Instr.distance      = Instr.distance + Instri.distance;
    Instr.centers       = Instr.centers + Instri.centers;
    Instr.gammaX        = Instr.gammaX + Instri.gammaX;
    Instr.gammaY        = Instr.gammaY + Instri.gammaY;
    Instr.detsizeHorz   = Instr.detsizeHorz + Instri.detsizeHorz;
    Instr.detsizeVert   = Instr.detsizeVert + Instri.detsizeVert;
    Instr.pixelsizeHorz	= Instr.pixelsizeHorz + Instri.pixelsizeHorz;
    Instr.pixelsizeVert	= Instr.pixelsizeVert + Instri.pixelsizeVert;
    Instr.numpixelsHorz	= Instr.numpixelsHorz + Instri.numpixelsHorz;
    Instr.numpixelsVert	= Instr.numpixelsVert + Instri.numpixelsVert;
    Instr.detpars       = Instr.detpars + Instri.detpars;
    Instr.imrotation    = 0;
    Instr.dettype       = '2a';
end
Instr.energy        = Instr.energy./length(XRDIMAGE.Calib.fnumber);
Instr.wavelength    = Instr.wavelength./length(XRDIMAGE.Calib.fnumber);
Instr.distance      = Instr.distance./length(XRDIMAGE.Calib.fnumber);
Instr.centers       = Instr.centers./length(XRDIMAGE.Calib.fnumber);
Instr.gammaX        = Instr.gammaX./length(XRDIMAGE.Calib.fnumber);
Instr.gammaY        = Instr.gammaY./length(XRDIMAGE.Calib.fnumber);
Instr.detsizeHorz   = Instr.detsizeHorz./length(XRDIMAGE.Calib.fnumber);
Instr.detsizeVert   = Instr.detsizeVert./length(XRDIMAGE.Calib.fnumber);
Instr.pixelsizeHorz	= Instr.pixelsizeHorz./length(XRDIMAGE.Calib.fnumber);
Instr.pixelsizeVert	= Instr.pixelsizeVert./length(XRDIMAGE.Calib.fnumber);
Instr.numpixelsHorz	= Instr.numpixelsHorz./length(XRDIMAGE.Calib.fnumber);
Instr.numpixelsVert	= Instr.numpixelsVert./length(XRDIMAGE.Calib.fnumber);
Instr.detpars       = Instr.detpars./length(XRDIMAGE.Calib.fnumber);

% XRDIMAGE.Image.samplename   = 'mishra_hr1_s1';
% FROOT   = {'mishra_hr1_s1_load11_pd'; ...
%     'mishra_hr1_s1_load12_pd'; ...
%     'mishra_hr1_s1_load13_pd'; ...
%     'mishra_hr1_s1_load14_pd'; ...
%     'mishra_hr1_s1_load15_pd'; ...
%     'mishra_hr1_s1_load16_pd' ; ...
%     'mishra_hr1_s1_load17_pd' ; ...
%     'mishra_hr1_s1_load18_pd' ; ...
%     'mishra_hr1_s1_load19_pd' ; ...
%     'mishra_hr1_s1_load20_pd' ; ...
%     'mishra_hr1_s1_load21_pd' ; ...
%     'mishra_hr1_s1_load22_pd' ; ...
%     'mishra_hr1_s1_load23_pd' ; ...
%     'mishra_hr1_s1_load24_pd' ; ...
%     'mishra_hr1_s1_load25_pd' ; ...
%     'mishra_hr1_s1_load26_pd' ; ...
%     'mishra_hr1_s1_load27_pd' ; ...
%     'mishra_hr1_s1_load28_pd' ; ...
%     'mishra_hr1_s1_load29_pd' ; ...
%     'mishra_hr1_s1_load30_pd' ; ...
%     'mishra_hr1_s1_load31_pd' ; ...
%     'mishra_hr1_s1_load32_pd' ; ...
%     };

% XRDIMAGE.Image.samplename   = 'mishra_hr2_s1';
% FROOT   = {'mishra_hr2_s1_load16_pd'; ...
%     'mishra_hr2_s1_load17_pd'; ...
%     'mishra_hr2_s1_load18_pd'; ...
%     'mishra_hr2_s1_load19_pd'; ...
%     'mishra_hr2_s1_load20_pd'; ...
%     'mishra_hr2_s1_load21_pd' ; ...
%     'mishra_hr2_s1_load22_pd' ; ...
%     'mishra_hr2_s1_load23_pd' ; ...
%     'mishra_hr2_s1_load24_pd' ; ...
%     'mishra_hr2_s1_load25_pd' ; ...
%     'mishra_hr2_s1_load26_pd' ; ...
%     'mishra_hr2_s1_load27_pd' ; ...
%     'mishra_hr2_s1_load28_pd' ; ...
%     'mishra_hr2_s1_load29_pd' ; ...
%     'mishra_hr2_s1_load30_pd' ; ...
%     'mishra_hr2_s1_load31_pd' ; ...
%     'mishra_hr2_s1_load32_pd' ; ...
%     'mishra_hr2_s1_load33_pd' ; ...
%     'mishra_hr2_s1_load34_pd' ; ...
%     'mishra_hr2_s1_load35_pd' ; ...
%     'mishra_hr2_s1_load36_pd' ; ...
%     'mishra_hr2_s1_load37_pd' ; ...
%     'mishra_hr2_s1_load38_pd' ; ...
%     'mishra_hr2_s1_load39_pd' ; ...
%     'mishra_hr2_s1_load40_pd' ; ...
%     'mishra_hr2_s1_load41_pd' ; ...
%     };

XRDIMAGE.Image.samplename   = 'mishra_hr3_s1';
FROOT   = {'mishra_hr3_s1_load16_pd'; ...
    'mishra_hr3_s1_load17_pd'; ...
    'mishra_hr3_s1_load18_pd'; ...
    'mishra_hr3_s1_load19_pd'; ...
    'mishra_hr3_s1_load20_pd'; ...
    'mishra_hr3_s1_load21_pd' ; ...
    'mishra_hr3_s1_load22_pd' ; ...
    'mishra_hr3_s1_load23_pd' ; ...
    'mishra_hr3_s1_load24_pd' ; ...
    'mishra_hr3_s1_load25_pd' ; ...
    'mishra_hr3_s1_load26_pd' ; ...
    'mishra_hr3_s1_load27_pd' ; ...
    'mishra_hr3_s1_load28_pd' ; ...
    'mishra_hr3_s1_load29_pd' ; ...
    'mishra_hr3_s1_load30_pd' ; ...
    'mishra_hr3_s1_load31_pd' ; ...
    'mishra_hr3_s1_load32_pd' ; ...
    'mishra_hr3_s1_load33_pd' ; ...
    'mishra_hr3_s1_load34_pd' ; ...
    'mishra_hr3_s1_load35_pd' ; ...
    'mishra_hr3_s1_load36_pd' ; ...
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATERIAL PARAMETERS - 
Material.num        = 1;
Material.lattparms  = 3.60;
Material.structure  = 'fcc';
Material.hkls       = load([Material.structure, '.hkls']);

%%% CALCULATE THEORETICAL TTH
[d, th] = PlaneSpacings(Material.lattparms, ...
    'cubic', Material.hkls', ...
    Instr.wavelength);
tth     = 2*th;
d_spacing_range = 0.03;
d_spacing_UB    = (1 + d_spacing_range)*d;
d_spacing_LB    = (1 - d_spacing_range)*d;

tth_UB  = 2.*asind(Instr.wavelength/2)./d_spacing_LB;
tth_LB  = 2.*asind(Instr.wavelength/2)./d_spacing_UB;

Material.tth        = tth;
Material.d_spacing  = d;
Material.pkidx      = {...
    [1] [2] [3]
    };

Material.numbounds  = length(Material.pkidx);
numpk   = 0;
for iii = 1:1:Material.numbounds
    Material.pkrange(:,iii)  = [ ...
        min(tth_LB(Material.pkidx{iii})); ...
        max(tth_UB(Material.pkidx{iii})); ...
        ];
    numpk   = numpk + length(Material.pkidx{iii});
end
Material.numpk  = numpk;
Material.pkbck  = 2;
Material.pkfunc = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PK FITTING OPTIONS
Analysis_Options.PkFuncOptions.pfunc_type	= 'pseudoVoigt';
Analysis_Options.PkFuncOptions.pbkg_order	= 2;
Analysis_Options.PkFitOptimizationOptions   = optimset(...
    'MaxIter', 5e5, ...
    'MaxFunEvals',3e5);

%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg    = 0;
Analysis_Options.save_polimg    = 0;
Analysis_Options.fits_spectra   = 1;
Analysis_Options.save_fits      = 1;
Analysis_Options.find_instrpars = 0;
Analysis_Options.save_instrpars = 0;
Analysis_Options.find_detpars	= 0;
Analysis_Options.generateESG    = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iiii = 1:1:length(FROOT)
    %%% INPUT PARAMETERS
    XRDIMAGE.Image.froot        = FROOT{iiii};     %%% FIND FNUM BASED ON THIS FILE NAME ROOT 32
    
    %%% GET IMAGE NUMBERS
    metadata_per_scan   = ReadSpecParFile(XRDIMAGE.Image.per_scan_metadata, 'Version', 'saxs_waxs_fmt_fastpar_v3');
    metadata_per_frame  = ReadSpecParFile(XRDIMAGE.Image.per_frame_metadata, 'Version' ,'mpe_ff_per_frame_v2');
    idx_froot           = ismember(metadata_per_scan.det3_fname, XRDIMAGE.Image.froot);
    
    XRDIMAGE.Image.fnumber      = metadata_per_scan.det3_fnum(idx_froot) - 1;
    
    switch XRDIMAGE.Image.corrected
        case 'sum'
            XRDIMAGE.Image.scan_nframes = 1;
        case {'cor', 'cor32'}
            XRDIMAGE.Image.scan_nframes = metadata_per_scan.scan_nframes(idx_froot);
            XRDIMAGE.Image.scan_mtr     = metadata_per_scan.scan_mtr(idx_froot);
            XRDIMAGE.Image.scan_ini     = metadata_per_scan.scan_ini(idx_froot);
            XRDIMAGE.Image.scan_fin     = metadata_per_scan.scan_fin(idx_froot);
        otherwise
            warning('num frames per scan could not be determined')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% GENERATE LIST OF FILES TO LOOK AT
    numimg  = length(XRDIMAGE.Image.fnumber);
    for iii = 1:1:numimg
        switch XRDIMAGE.Image.corrected
            case 'sum'
                fname   = sprintf('%s_%06d.%s.%s', ...
                    XRDIMAGE.Image.froot, XRDIMAGE.Image.fnumber(iii), XRDIMAGE.Image.fext, XRDIMAGE.Image.corrected);
                pfname{iii, 1}  = fullfile(XRDIMAGE.Image.pname, ...
                    XRDIMAGE.Image.samplename, XRDIMAGE.Image.fext, fname);
            case {'cor', 'cor32'}
                for jjj = 1:1:XRDIMAGE.Image.scan_nframes
                    fname   = sprintf('%s_%06d.%s_frame_%d.%s', ...
                        XRDIMAGE.Image.froot, XRDIMAGE.Image.fnumber(iii), XRDIMAGE.Image.fext, jjj, XRDIMAGE.Image.corrected);
                    pfname{iii, jjj}    = fullfile(XRDIMAGE.Image.pname, ...
                        XRDIMAGE.Image.samplename, XRDIMAGE.Image.fext, fname);
                end
            otherwise
                warning('check XRDIMAGE.Image.corrected paramter')
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LOAD XRD IMAGES & GENERATE POLIMG FILE
    if Analysis_Options.fits_spectra
        for iii = 1:1:numimg
            scan_mtr    = XRDIMAGE.Image.scan_mtr{iii};
            switch scan_mtr
                case 'aero'
                    omega_grid  = linspace(XRDIMAGE.Image.scan_ini(iii), ...
                        XRDIMAGE.Image.scan_fin(iii), ...
                        XRDIMAGE.Image.scan_nframes(iii)+1);
                    omega_grid  = (omega_grid(1:end-1) + omega_grid(2:end))./2;
                    chi_grid    = 0.*omega_grid;
                case 'prrot'
                    warning('we have nothing here; maybe implement a flag that says prrot is angle or not?')
            end
            
            for jjj = 1:1:XRDIMAGE.Image.scan_nframes(iii)
                tic
                disp('###########################')
                disp(sprintf('Looking at %s', pfname{iii, jjj}))
                disp('###########################')
                
                pfname_polimg   = sprintf('%s.polimg.mat', pfname{iii,jjj});
                pfname_pkfit    = sprintf('%s.pkfit.mat', pfname{iii,jjj});
                
                polimg  = load(pfname_polimg);
                
                CakePrms    = polimg.cakeprms;
                Instr       = polimg.instr;
                polimg      = polimg.polimg;
                
                figure(2)
                subplot(1,2,1)
                imagesc(log(abs(polimg.intensity_in_tth_grid))), axis square tight
                hold off
                
                subplot(1,2,2)
                plot(polimg.tth_grid, polimg.intensity_in_tth_grid)
                hold off
                
                disp('###########################')
                disp(sprintf('Fitting peaks in %s', pfname{iii,1}))
                disp('###########################')
                
                for kkk = 1:1:CakePrms.bins(1)
                    fprintf('Looking at azimuthal bin %d of %d\n', kkk, CakePrms.bins(1))
                    
                    x   = polimg.tth_grid;
                    y   = polimg.intensity_in_tth_grid(kkk,:);
                    
                    figure(11)
                    subplot(1,2,1)
                    plot(x, y, 'k.')
                    hold on
                    plot(tth, mean(y), 'g^')
                    axis([min(x) max(x) min(y) max(y)+100])
                    xlabel('radial distance (mm)')
                    ylabel('intensity (arb. units)')
                    title(['bin number : ', num2str(kkk)])
                    
                    for mmm = 1:1:Material.numbounds
                        numpks  = length(Material.pkidx{mmm});
                        fprintf('Looking at bound number %d of %d with %d peaks\n', mmm, Material.numbounds, numpks);
                        
                        pkrange = Material.pkrange(:,mmm);
                        
                        idx = find(x >= pkrange(1) & x <= pkrange(2));
                        xr  = x(idx)';
                        yr  = y(idx)';
                        
                        if kkk == 1
                            idx_max = [];
                            peakdet_thresh  = 0.5;
                            while (length(idx_max) ~= numpks) && (peakdet_thresh > 0)
                                [idx_max, idx_min]  = peakdet(yr, peakdet_thresh, xr);
                                peakdet_thresh      = peakdet_thresh - 0.1;
                            end
                            
                            %                             pr0 = [];
                            %                             if length(idx_max) == numpks
                            %                                 for nnn = 1:1:numpks
                            %                                 % pr0 =
                            %                                 end
                            %                             else
                            %                                 for nnn = 1:1:numpks
                            %                                     % pr0 =
                            %                                 end
                            %                             end
                            
                            pkidx   = Material.pkidx{mmm};
                            pr0     = [];
                            pr0_LB  = [];
                            pr0_UB  = [];
                            for nnn = 1:1:numpks
                                pr0 = [ ...
                                    pr0;
                                    idx_max(nnn,2)/10; 
                                    0.05; 
                                    0.5; 
                                    tth(pkidx(nnn)); 
                                    ];
                                
                                pr0_LB  = [ ...
                                    pr0_LB; ...
                                    0; ...
                                    0; ...
                                    0; ...
                                    tth(pkidx(nnn))-0.5; ...
                                    ];
                                
                                pr0_UB  = [ ...
                                    pr0_UB; ...
                                    inf; ...
                                    inf; ...
                                    1; ...
                                    tth(pkidx(nnn))+0.5; ...
                                    ];
                            end
                            
                            pr0 = [pr0; ...
                                0; ...
                                yr(1); ...
                                ];
                            
                            pr0_LB  = [pr0_LB; ...
                                -inf; ...
                                -inf; ...
                                ];
                            
                            pr0_UB  = [pr0_UB; ...
                                inf; ...
                                inf; ...
                                ];
                        else
                            pr0     = pr_previous_azimuth{mmm};
                            pr0_LB  = pr_LB_previous_azimuth{mmm};
                            pr0_UB  = pr_UB_previous_azimuth{mmm};
                        end
                        
                        pkpars.pfunc_type   = Analysis_Options.PkFuncOptions.pfunc_type;
                        pkpars.pbkg_order   = Analysis_Options.PkFuncOptions.pbkg_order;
                        pkpars.xdata        = xr;
                        
                        [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc_switch, pr0, pkpars, yr, ...
                            pr0_LB, pr0_UB, Analysis_Options.PkFitOptimizationOptions);
                        y0	= pfunc_switch(pr0, pkpars);
                        yf	= pfunc_switch(pr, pkpars);
                        
                        % [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc, pr0, xr, yr, ...
                        %     [], [], Analysis_Options.PkFitOptions);
                        % y0  = pfunc(pr0,xr);
                        % yf  = pfunc(pr,xr);
                        
                        %%% UPDATE PREVIOUS AZIMUTH PR FOR NEXT AZIM FITING
                        pr_previous_azimuth{mmm}    = pr;
                        pr_LB_previous_azimuth{mmm} = pr0_LB;
                        pr_UB_previous_azimuth{mmm} = pr0_UB;
                        
                        figure(11)
                        subplot(1,2,1)
                        plot(xr, yr, 'b.')
                        plot(xr, y0, 'r-')
                        plot(xr, yf, 'g-')
                        
                        subplot(1,2,2)
                        plot(xr, yr, 'b.')
                        hold on
                        plot(xr, y0, 'r-')
                        plot(xr, yf, 'g-')
                        xlabel('radial distance (mm)')
                        ylabel('intensity (arb. units)')
                        title(['bound number : ', num2str(mmm)])
                        hold off
                        
                        %%% MAPPING NEEDS TO BE UPDATED WITH A SWITCH
                        % pro  = pkfitResultMapping(pkpars, pr);
                        
                        pkfit.amp(kkk,mmm)  = pr(1);
                        pkfit.fwhm(kkk,mmm) = pr(2);
                        pkfit.mix(kkk,mmm)  = pr(3);
                        pkfit.rho(kkk,mmm)  = pr(4);
                        pkfit.bkg{kkk,mmm}  = pr(5:end);
                        pkfit.rsn(kkk,mmm)  = rsn;
                        pkfit.ef(kkk,mmm)   = ef;
                        pkfit.rwp(kkk,mmm)  = ErrorRwp(yr, yf);
                    end
                    figure(11)
                    subplot(1,2,1)
                    hold off
                end
                
                if Analysis_Options.save_fits
                    disp('###########################')
                    disp(sprintf('Saving peak fits in %s\n', pfname_pkfit))
                    save(pfname_pkfit, 'pkfit')
                else
                    disp('###########################')
                    disp(sprintf('Not saving peak fits for %s\n', pfname{iii, jjj}))
                end
                toc
            end
        end
    end
end