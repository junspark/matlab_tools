clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS FROM CALIBRANT 
DataReductionPrms.Dsam      = 1874.565;
DataReductionPrms.Energy    = 65.351;
DataReductionPrms.Lambda    = keV2Angstrom(DataReductionPrms.Energy);
DataReductionPrms.x0        = 2279.502;
DataReductionPrms.y0        = 2118.644;
DataReductionPrms.TiltPlane = 55.933;
DataReductionPrms.InPlane   = 1.451;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CAKING PARAMETERS
DataReductionPrms.PixSize   = 0.2;
DataReductionPrms.ETAStart  = -170;
DataReductionPrms.ETAEnd    = -114;
DataReductionPrms.RHOInner  = 500;
DataReductionPrms.RHOOuter  = 2372;
DataReductionPrms.ETABins   = 13;
DataReductionPrms.RHOBins   = 2048;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALIBRANT GE FILE NAME
pname   = 'X:\balogh_march14\calibration';
fname   = 'Ceo2_calibr1_00030.ge4';
pfname  = fullfile(pname, fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATERIAL INFORMATION
latticeParms    = 5.411651;
hkls        = load('fcc.hkls');
[d0, th0]   = PlaneSpacings(latticeParms, 'cubic', hkls', DataReductionPrms.Lambda);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tth0    = 2.*th0;
x0  = DataReductionPrms.Dsam.*tand(tth0);

deta    = (DataReductionPrms.ETAEnd - DataReductionPrms.ETAStart)/DataReductionPrms.ETABins;
eta     = DataReductionPrms.ETAStart:deta:(DataReductionPrms.ETAEnd-deta);
dr      = (DataReductionPrms.RHOOuter - DataReductionPrms.RHOInner)/(DataReductionPrms.RHOBins - 1);

%%% LOAD SPR FILE
SPRName	= fullfile(pname, ...
    [fname, '.spr']);

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
imagesc(spr_data);
title(fname, 'Interpreter', 'none')
xticklabel  = (str2num(get(gca, 'XTickLabel'))./dr + DataReductionPrms.RHOInner).*DataReductionPrms.PixSize;
xticklabel  = atand(xticklabel./DataReductionPrms.Dsam);
xticklabel  = num2str(xticklabel, '%10.2f');
yticklabel  = (str2num(get(gca, 'YTickLabel')).*deta + DataReductionPrms.ETAStart);
yticklabel  = num2str(yticklabel, '%10.2f');
set(gca, 'XTickLabel', xticklabel)
set(gca, 'YTickLabel', yticklabel)
xlabel('radial position - 2\theta (deg)')
ylabel('azimuthal position \eta (deg)')

%%% FIT PEAKS
x   = linspace(DataReductionPrms.RHOInner, DataReductionPrms.RHOOuter, DataReductionPrms.RHOBins);
x   = x.*DataReductionPrms.PixSize;
for j = 1:1:DataReductionPrms.ETABins
    y   = spr_data(j,:);
    
    figure(1)
    clf
    plot(x,y, 'k.')
    hold on
    plot(x0, ones(length(x0),1), 'r^')
    
    for k = 1:1:17%length(tth0)
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
        title([fname, '-azimuth: ', num2str(eta(j)), '-peak number:', num2str(k)], 'Interpreter', 'none')
        pause(0.1)
        
        tth(j,k) = atand(pjk(4)/DataReductionPrms.Dsam);
    end
end
d           = 0.5*DataReductionPrms.Lambda./sind(tth./2);
psstrain    = (d - repmat(d0(1:17), 13, 1))./repmat(d0(1:17), 13, 1);

figure, 
imagesc(psstrain)
colorbar vert
return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% USER INPUT
% %%% LOCATION AND NAME OF DARKFILE
% DataReductionPrms.DarkPath  = '.';
% DataReductionPrms.DarkName  = { ...
%     'CeO2_00680.ge2'; ...
%     };
% 
% %%% LOCATION AND NAME OF A HYDRA IMAGE
% DataReductionPrms.ImagePath     = '.';
% DataReductionPrms.ImageNames    = { ...
%     'Ceo2_calibr1_00030.ge1'; ...
%     };
% 
% %%% LOCATION AND NAME OF FIT2D MACRO
% % DataReductionPrms.MacroPath = '.';
% % DataReductionPrms.MacroName = 'fit2dmacro';
% % DataReductionPrms.SPRPath   = '.';
% 
% %%% LOCATION OF TTH RESULTS
% DataReductionPrms.TTHPath   = '.';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% deta    = (DataReductionPrms.ETAEnd - DataReductionPrms.ETAStart)/DataReductionPrms.ETABins;
% eta     = DataReductionPrms.ETAStart:deta:(DataReductionPrms.ETAEnd-deta);
% 
% tth0    = [ ...
%     4.07; ...
%     5.76; ...
%     7.06; ...
%     ];
% x0  = DataReductionPrms.Dsam.*tand(tth0);
% 
% %%% GENERATE SPR FILES
% [MacroName, LogFileName] = WriteFit2DCakeMacro(DataReductionPrms);
% disp('running fit2d_12_077_i686_WXP')
% disp(['run macro ', MacroName])
% ! C:\Users\parkjs\Applications\fit2d_12_077_i686_WXP.exe
% 
% return
% %%% NOW FIT
% ni  = length(DataReductionPrms.ImageNames);
% for i = 1:1:ni
%     %%% LOAD SPR FILE
%     SPRName	= fullfile(DataReductionPrms.SPRPath, ...
%         [DataReductionPrms.ImageNames{i}, '.spr']);
%     
%     spr_data    = [];
%     fid_spr = fopen(SPRName, 'r');
%     while ~feof(fid_spr)
%         lindata = fgetl(fid_spr);
%         if isempty(strfind(lindata, 'Start pixel'))
%             spr_data    = [spr_data; str2num(lindata)];
%         end
%     end
%     fclose(fid_spr);
%     figure(100)
%     imagesc(spr_data)
%     title(DataReductionPrms.ImageNames{i}, 'Interpreter', 'none')
%     
%     %%% FIT PEAKS
%     x   = linspace(DataReductionPrms.RHOInner, DataReductionPrms.RHOOuter, DataReductionPrms.RHOBins);
%     x   = x.*DataReductionPrms.PixSize;
%     for j = 1:1:DataReductionPrms.ETABins
%         y   = spr_data(j,:);
%         
%         figure(1)
%         clf
%         plot(x,y, 'k.')
%         hold on
%         plot(x0, ones(3,1), 'r^')
%         
%         for k = 1:1:3
%             idx = (x < (x0(k) + 2)) & (x > (x0(k) - 2));
%             xjk = x(idx)';
%             yjk = y(idx)';
%             
%             %%% SIMPLEST FIT
%             pjk0    = [max(yjk)/2; 0.5; 0.025; x0(k); 0; 0];
%             
%             pjk0    = [pjk0; 0; 30];
%             yjk0    = pfunc(pjk0, xjk);
%             
%             pjk     = lsqcurvefit(@pfunc, pjk0, xjk, yjk);
%             yjkf    = pfunc(pjk, xjk);
% 
%             plot(xjk, yjk, 'b.')
%             plot(xjk, yjk0, 'r-')
%             plot(xjk, yjkf, 'g-')
%             xlabel('radial position (mm)')
%             ylabel('intensity (arb. units)')
%             title([DataReductionPrms.ImageNames{i}, '-azimuth: ', num2str(eta(j)), '-peak number:', num2str(k)], 'Interpreter', 'none')
%             pause(0.1)
%             
%             tth(j,k) = atand(pjk(4)/DataReductionPrms.Dsam);
%         end
%     end
%     
%     fname_tth   = [DataReductionPrms.ImageNames{i}, '.tth'];
%     pfname_tth  = fullfile(DataReductionPrms.TTHPath, fname_tth);
%     
%     fid_tth = fopen(pfname_tth, 'w');
%     fprintf(fid_tth, '%% azimuth (deg) tth_110 (deg) tth_200 (deg) tth_211 (deg)');
%     fprintf(fid_tth, '%f\t%f\t%f\t%f\n', [eta' tth]');
%     fclose(fid_tth);
% end