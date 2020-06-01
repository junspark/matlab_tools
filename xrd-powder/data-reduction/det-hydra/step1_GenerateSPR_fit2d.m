clear all
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT
%%% LOCATION AND NAME OF STRAIN IMAGES - PROVIDE A LIST OF IMAGE FILE NAMES
DataReductionPrms.ImagePath     = 'O:\motta_nov13\APS2_corrected';
fnum    = 12:1:6481;
fstem   = 'APS2_';
fext    = 'ge3.sum';
for i = 1:1:length(fnum)
    DataReductionPrms.ImageNames{i,1}   = sprintf([fstem, '%05d.%s'], fnum(i), fext);
end

%%% LOCATION AND NAME OF FIT2D MACRO
DataReductionPrms.MacroPath = '.';
DataReductionPrms.MacroName = 'fit2dmacro.ge3';
DataReductionPrms.SPRPath   = '.';

%%% LOCATION OF TTH RESULTS
DataReductionPrms.TTHPath   = '.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS FROM CALIBRANT 
DataReductionPrms.Dsam      = 2150.805;
DataReductionPrms.Energy    = 51.996;
DataReductionPrms.Lambda    = keV2Angstrom(DataReductionPrms.Energy);
DataReductionPrms.x0        = 2264.611;
DataReductionPrms.y0        = 2134.987;
DataReductionPrms.TiltPlane = -130.169;
DataReductionPrms.InPlane   = 1.329;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CAKING PARAMETERS
DataReductionPrms.PixSize   = 0.2;
DataReductionPrms.ETAStart  = -150;
DataReductionPrms.ETAEnd    = -120;
DataReductionPrms.RHOInner  = 650;
DataReductionPrms.RHOOuter  = 2600;
DataReductionPrms.ETABins   = 6;
DataReductionPrms.RHOBins   = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deta    = (DataReductionPrms.ETAEnd - DataReductionPrms.ETAStart)/DataReductionPrms.ETABins;
eta     = DataReductionPrms.ETAStart:deta:(DataReductionPrms.ETAEnd-deta);

hkls    = load('fcc.hkls')';
[d, th] = PlaneSpacings(5.411651, 'cubic', hkls, DataReductionPrms.Lambda);

tth0    = 2.*th(1:10);
x0      = DataReductionPrms.Dsam.*tand(tth0);

%%% GENERATE SPR FILES
[MacroName, LogFileName] = WriteFit2DCakeMacro2(DataReductionPrms);
disp('running fit2d_12_077_i686_WXP')
disp(['run macro ', MacroName])
! C:\Users\parkjs\Applications\fit2d_12_077_i686_WXP.exe

return
%%% NOW FIT
ni  = length(DataReductionPrms.ImageNames);
for i = 1:1:ni
    %%% LOAD SPR FILE
    SPRName	= fullfile(DataReductionPrms.SPRPath, ...
        [DataReductionPrms.ImageNames{i}, '.spr']);
    
    spr_data    = [];
    fid_spr = fopen(SPRName, 'r');
    while ~feof(fid_spr)
        lindata = fgetl(fid_spr);
        if isempty(strfind(lindata, 'Start pixel'))
            spr_data    = [spr_data; str2num(lindata)];
        end
    end
    fclose(fid_spr);
    figure(100)
    imagesc(spr_data)
    title(DataReductionPrms.ImageNames{i}, 'Interpreter', 'none')
    pause(0.5)
    %%% FIT PEAKS
    x   = linspace(DataReductionPrms.RHOInner, DataReductionPrms.RHOOuter, DataReductionPrms.RHOBins);
    x   = x.*DataReductionPrms.PixSize;
    for j = 1:1:DataReductionPrms.ETABins
        y   = spr_data(j,:);
        
        figure(1)
        clf
        plot(x,y, 'k.')
        hold on
        plot(x0, ones(length(tth0),1), 'r^')
        
        for k = 1:1:length(tth0)
            idx = (x < (x0(k) + 10)) & (x > (x0(k) - 10));
            xjk = x(idx)';
            yjk = y(idx)';
            
            %%% SIMPLEST FIT
            pjk0    = [max(yjk)/2; 0.5; 0.025; x0(k); 0; 0];
            
            pjk0    = [pjk0; 0; 30];
            yjk0    = pfunc(pjk0, xjk);
            
            pjk     = lsqcurvefit(@pfunc, pjk0, xjk, yjk);
            yjkf    = pfunc(pjk, xjk);

            plot(xjk, yjk, 'b.')
            plot(xjk, yjk0, 'r-')
            plot(xjk, yjkf, 'g-')
            xlabel('radial position (mm)')
            ylabel('intensity (arb. units)')
            title([DataReductionPrms.ImageNames{i}, '-azimuth: ', num2str(eta(j)), '-peak number:', num2str(k)], 'Interpreter', 'none')
            pause(0.1)
            
            tth(j,k) = atand(pjk(4)/DataReductionPrms.Dsam);
        end
    end
    
    fname_tth   = [DataReductionPrms.ImageNames{i}, '.tth'];
    pfname_tth  = fullfile(DataReductionPrms.TTHPath, fname_tth);
    
    fid_tth = fopen(pfname_tth, 'w');
    fprintf(fid_tth, '%% azimuth (deg) tth_n (deg)');
    fprintf(fid_tth, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', [eta' tth]');
    fclose(fid_tth);
end
pseudostrain    = sind(repmat(tth0, DataReductionPrms.ETABins, 1)./2)./sind(tth./2) - 1;
figure,
imagesc(pseudostrain)
colorbar vert