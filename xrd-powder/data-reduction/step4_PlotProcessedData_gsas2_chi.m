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

data    = load(fname_reprocessed);
numpts  = length(data.amp);
numpks  = data.numpks;

hkls    = load('fcc.hkls'); %%% NOTE THIS PEAK ID

%%% PLOT VARIOUS MEASURED QUANTITIES FOR PARTICULAR DIRECTIONS
%%% THESE ARE ALL NOMINAL LD AND TD
%%% ETA 0 / 180 == TD; ETA 90 / 270 == LD
[~, idx_TD1]    = sort(abs(data.eta{1} - 0));
[~, idx_TD2]    = sort(abs(data.eta{1} - 180));
[~, idx_LD1]    = sort(abs(data.eta{1} - 90));
[~, idx_LD2]    = sort(abs(data.eta{1} - 270));

for iii = 1:1:numpts
    for jjj = 1:1:numpks
        %%%
        eTD{jjj}(iii,1)  = data.lattice_strain{iii}(idx_TD1(1),jjj);
        eTD{jjj}(iii,2)  = data.lattice_strain{iii}(idx_TD2(2),jjj);
        
        eLD{jjj}(iii,1)  = data.lattice_strain{iii}(idx_LD1(1),jjj);
        eLD{jjj}(iii,2)  = data.lattice_strain{iii}(idx_LD2(2),jjj);
        eAllAzim{jjj}(iii,:)    = data.lattice_strain{iii}(:,jjj);
        
        %%%
        aTD{jjj}(iii,1)  = data.amp{iii}(idx_TD1(1),jjj);
        aTD{jjj}(iii,2)  = data.amp{iii}(idx_TD2(2),jjj);
        
        aLD{jjj}(iii,1)  = data.amp{iii}(idx_LD1(1),jjj);
        aLD{jjj}(iii,2)  = data.amp{iii}(idx_LD2(2),jjj);
        aAllAzim{jjj}(iii,:)    = data.amp{iii}(:,jjj);
        
        %%%
        fTD{jjj}(iii,1)  = data.fwhm{iii}(idx_TD1(1),jjj);
        fTD{jjj}(iii,2)  = data.fwhm{iii}(idx_TD2(2),jjj);
        
        fLD{jjj}(iii,1)  = data.fwhm{iii}(idx_LD1(1),jjj);
        fLD{jjj}(iii,2)  = data.fwhm{iii}(idx_LD2(2),jjj);
        fAllAzim{jjj}(iii,:)    = data.fwhm{iii}(:,jjj);
    end
end

%%%% STRAIN
xplot   = XRDIMAGE.Image.fnumber;          %% THIS IS LOAD IN N
for iii = 1:1:numpks
    figure(iii)
    hold off
    plot(xplot, eTD{iii}(:,1), 'r-')
    hold on
    plot(xplot, eTD{iii}(:,2), 'r--')
    plot(xplot, eLD{iii}(:,1), 'b-')
    plot(xplot, eLD{iii}(:,2), 'b--')
    xlabel('image number (-)')
    ylabel('lattice strain (-)')
    hold off
    % axis([-100 1300 -2e-3 8e-3])
    title(sprintf('hkl = %d%d%d', hkls(iii,1), hkls(iii,2), hkls(iii,3)))
    grid on
    
    fname_fig   = sprintf('%s.strain_load.hkl.%d%d%d.png', ...
        XRDIMAGE.Image.fbase, hkls(iii,1), hkls(iii,2), hkls(iii,3));
    saveas(gcf, fname_fig, 'png')
end
close all

% xplot   = data.enc(4,:);          %% THIS IS DISP IN mm
% for iii = 1:1:numpks
%     
%     figure(iii)
%     hold off
%     imagesc(eAllAzim{iii}')
%     colorbar vert
%     caxis([-2e-3 8e-3])
%     colormap jet
%     axis square tight
%     xlabel('disp (mm)')
%     ylabel('azimuth (deg)')
%     title(sprintf('hkl = %d%d%d', hkls(iii,1), hkls(iii,2), hkls(iii,3)))
%     
%     aaa = gca;
%     idx = get(aaa, 'YTick');
%     set(gca, 'YTickLabel', num2str(floor(data.eta{iii}(idx))))
%     
%     idx = get(aaa, 'XTick');
%     set(gca, 'XTickLabel', num2str(round(xplot(idx)'.*1000)./1000))    
%     
%     fname_fig   = sprintf('%s.azim_disp_strain.hkl.%d%d%d.png', ...
%         XRDIMAGE.Image.fbase, hkls(iii,1), hkls(iii,2), hkls(iii,3));
%     saveas(gcf, fname_fig, 'png')
% end
% close all
% 
% %%%% FWHM
% xplot   = data.ev(2,:);          %% THIS IS LOAD IN N
% for iii = 1:1:numpks
%     figure(iii+100)
%     hold off
%     plot(xplot, fTD{iii}(:,1), 'r-')
%     hold on
%     plot(xplot, fTD{iii}(:,2), 'r--')
%     plot(xplot, fLD{iii}(:,1), 'b-')
%     plot(xplot, fLD{iii}(:,2), 'b--')
%     xlabel('load (N)')
%     ylabel('fwhm')
%     hold off
%     axis([-100 1300 0 1])
%     title(sprintf('hkl = %d%d%d', hkls(iii,1), hkls(iii,2), hkls(iii,3)))
%     grid on
%     
%     fname_fig   = sprintf('%s.fwhm_load.hkl.%d%d%d.png', ...
%         XRDIMAGE.Image.fbase, hkls(iii,1), hkls(iii,2), hkls(iii,3));
%     saveas(gcf, fname_fig, 'png')
% end
% close all
% 
% xplot   = data.enc(4,:);          %% THIS IS DISP IN mm
% for iii = 1:1:numpks
%     
%     figure(iii)
%     hold off
%     imagesc(fAllAzim{iii}')
%     colorbar vert
%     caxis([0 1.0])
%     colormap jet
%     axis square tight
%     xlabel('disp (mm)')
%     ylabel('azimuth (deg)')
%     title(sprintf('hkl = %d%d%d', hkls(iii,1), hkls(iii,2), hkls(iii,3)))
%     
%     aaa = gca;
%     idx = get(aaa, 'YTick');
%     set(gca, 'YTickLabel', num2str(floor(data.eta{iii}(idx))))
%     
%     idx = get(aaa, 'XTick');
%     set(gca, 'XTickLabel', num2str(round(xplot(idx)'.*1000)./1000))    
%     
%     fname_fig   = sprintf('%s.azim_disp_fwhm.hkl.%d%d%d.png', ...
%         XRDIMAGE.Image.fbase, hkls(iii,1), hkls(iii,2), hkls(iii,3));
%     saveas(gcf, fname_fig, 'png')
% end
% close all
% 
% %%% AMPLITUDE
% xplot   = data.enc(4,:);          %% THIS IS DISP IN mm
% for iii = 1:1:numpks
%     
%     %%% APPROXIMATE INTGRATED PEAK INTENSITY
%     fff = log(aAllAzim{iii}.*fAllAzim{iii});
%     
%     figure(iii)
%     hold off
%     imagesc(fff')
%     colorbar vert
%     caxis([0 5])
%     colormap jet
%     axis square tight
%     xlabel('disp (mm)')
%     ylabel('azimuth (deg)')
%     title(sprintf('hkl = %d%d%d', hkls(iii,1), hkls(iii,2), hkls(iii,3)))
%     
%     aaa = gca;
%     idx = get(aaa, 'YTick');
%     set(gca, 'YTickLabel', num2str(floor(data.eta{iii}(idx))))
%     
%     idx = get(aaa, 'XTick');
%     set(gca, 'XTickLabel', num2str(round(xplot(idx)'.*1000)./1000))    
%     
%     fname_fig   = sprintf('%s.azim_disp_amp.hkl.%d%d%d.png', ...
%         XRDIMAGE.Image.fbase, hkls(iii,1), hkls(iii,2), hkls(iii,3));
%     saveas(gcf, fname_fig, 'png')
% end
% close all