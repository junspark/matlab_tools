clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT
%%% LOCATION AND NAME OF DARKFILE
DataReductionPrms.DarkPath  = 'C:\Users\parkjs\Documents\GitHub\matlab_tools\xrd-powder\data-reduction\example\APS';
DataReductionPrms.DarkName  = { ...
    'LaB6_04115.ge2'; ...
    };

%%% LOCATION AND NAME OF STRAIN IMAGES - PROVIDE A LIST OF IMAGE FILE NAMES
DataReductionPrms.ImagePath     = 'C:\Users\parkjs\Documents\GitHub\matlab_tools\xrd-powder\data-reduction\example\APS';
DataReductionPrms.ImageNames    = { ...
    'LaB6_04116.ge2'; ...
    'LaB6_04117.ge2'; ...
    };

%%% LOCATION AND NAME OF FIT2D MACRO
DataReductionPrms.MacroPath = '.';
DataReductionPrms.MacroName = 'fit2dmacro';
DataReductionPrms.SPRPath   = 'C:\Users\parkjs\Documents\GitHub\matlab_tools\xrd-powder\data-reduction\example\APS';

%%% LOCATION OF TTH RESULTS
DataReductionPrms.TTHPath   = '.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS FROM CALIBRANT 
DataReductionPrms.Dsam      = 1778.581;
DataReductionPrms.Energy    = 86;
DataReductionPrms.Lambda    = keV2Angstrom(DataReductionPrms.Energy);
DataReductionPrms.x0        = 1029.015;
DataReductionPrms.y0        = 1026.167;
DataReductionPrms.TiltPlane = -21.842;
DataReductionPrms.InPlane   = 0.271;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CAKING PARAMETERS
DataReductionPrms.PixSize   = 0.2;
DataReductionPrms.ETAStart  = 0;
DataReductionPrms.ETAEnd    = 360;
DataReductionPrms.RHOInner  = 200;
DataReductionPrms.RHOOuter  = 1000;
DataReductionPrms.ETABins   = 72;
DataReductionPrms.RHOBins   = 1500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deta    = (DataReductionPrms.ETAEnd - DataReductionPrms.ETAStart)/DataReductionPrms.ETABins;
eta     = DataReductionPrms.ETAStart:deta:(DataReductionPrms.ETAEnd-deta);

tth0    = [ ...
    4.07; ...
    5.76; ...
    7.06; ...
    ];
x0  = DataReductionPrms.Dsam.*tand(tth0);

%%% GENERATE SPR FILES
[MacroName, LogFileName] = WriteFit2DCakeMacro(DataReductionPrms);
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
    
    %%% FIT PEAKS
    x   = linspace(DataReductionPrms.RHOInner, DataReductionPrms.RHOOuter, DataReductionPrms.RHOBins);
    x   = x.*DataReductionPrms.PixSize;
    for j = 1:1:DataReductionPrms.ETABins
        y   = spr_data(j,:);
        
        figure(1)
        clf
        plot(x,y, 'k.')
        hold on
        plot(x0, ones(3,1), 'r^')
        
        for k = 1:1:3
            idx = (x < (x0(k) + 2)) & (x > (x0(k) - 2));
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
    fprintf(fid_tth, '%% azimuth (deg) tth_110 (deg) tth_200 (deg) tth_211 (deg)');
    fprintf(fid_tth, '%f\t%f\t%f\t%f\n', [eta' tth]');
    fclose(fid_tth);
end