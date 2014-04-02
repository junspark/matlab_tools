clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS FROM CALIBRANT 
DataReductionPrms.Dsam      = 1875.359;
DataReductionPrms.Energy    = 65.351;
DataReductionPrms.Lambda    = keV2Angstrom(DataReductionPrms.Energy);
DataReductionPrms.x0        = 2263.547;
DataReductionPrms.y0        = 2138.371;
DataReductionPrms.TiltPlane = 78.79421;
DataReductionPrms.InPlane   = -0.190839;
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
%%% GE FILE NAME
DataReductionPrms.ImagePath = 'O:\balogh_march14\GE1';
DataReductionPrms.SPRPath   = 'W:\balogh_march14\GE1';
fstem       = 'ss316block3b';
imnum_ini   = 155;
imnum_fin   = 158;
ct  = 1;
for i = imnum_ini:1:imnum_fin
    fname   = sprintf('%s_%05d.ge1.sum', fstem, i);
    pfname  = fullfile(DataReductionPrms.ImagePath, fname);
    if exist(pfname, 'file') == 2
        DataReductionPrms.ImageNames{ct,1}   = fname;
    end
    ct  = ct + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOCATION AND NAME OF FIT2D MACRO
DataReductionPrms.MacroPath = '.';
DataReductionPrms.MacroName = 'fit2dmacro.ge1';
[MacroName, LogFileName] = WriteFit2DCakeMacro2(DataReductionPrms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! C:\Users\parkjs\Applications\fit2d_12_077_i686_WXP.exe