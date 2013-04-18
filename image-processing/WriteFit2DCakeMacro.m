function [pfMacro, pfLog] = WriteFit2DCakeMacro(DataReductionPrms)
% Creates fit2d macro file
DarkPath    = DataReductionPrms.DarkPath;
DarkName    = DataReductionPrms.DarkName;
DarkName    = fullfile(DarkPath, DarkName{1});

ImagePath   = DataReductionPrms.ImagePath;
ImageNames  = DataReductionPrms.ImageNames;

MacroPath   = DataReductionPrms.MacroPath;
MacroName   = DataReductionPrms.MacroName;

SPRPath     = DataReductionPrms.SPRPath;

Dsam        = DataReductionPrms.Dsam;
Lambda      = DataReductionPrms.Lambda;
centerX     = DataReductionPrms.x0;
centerY     = DataReductionPrms.y0;
TiltPlane   = DataReductionPrms.TiltPlane;
InPlane     = DataReductionPrms.InPlane;

pixSize     = DataReductionPrms.PixSize*1e3;

etaStart    = DataReductionPrms.ETAStart;
etaEnd      = DataReductionPrms.ETAEnd;
rhoStart    = DataReductionPrms.RHOInner;
rhoEnd      = DataReductionPrms.RHOOuter;
numBinsA    = DataReductionPrms.ETABins;
numBinsR    = DataReductionPrms.RHOBins;

ni	= length(ImageNames);
SPRNames    = cell(ni,1);
for i = 1:1:ni
    SPRNames{i}     = fullfile(SPRPath, [ImageNames{i}, '.spr']);
    ImageNames{i}   = fullfile(ImagePath, ImageNames{i});
end

% create macro file to write
pfMacro = fullfile(MacroPath, MacroName);
pfMacro = [pfMacro, '.mac'];
pfLog   = fullfile(MacroPath, MacroName);
pfLog   = [pfLog, '.f2dlog'];

fid = fopen(pfMacro, 'w');

% header information
header = {...
    '%!*\ BEGINNING OF GUI MACRO FILE'
    '%!*\'
    '%!*\ This is a comment line'};
for i = 1:1:ni
    header{end+1}   = ['%!*\ fit2d macro for ', ImageNames{i}];
end
header{end+1}   = '%!*\';
header{end+1}   = 'OPEN LOG FILE';
header{end+1}   = 'NO';
header{end+1}   = pfLog;
header{end+1}   = 'EXIT';
header{end+1}   = 'POWDER DIFFRACTION (2-D)';

% write header of the macro file
for i = 1:length(header)
    fprintf(fid, '%s\n', header{i});
end

% write body of the macro file
body = [];
for i = 1:1:ni
    if i == 1
        body{i} = {...
            'INPUT'
            [ImageNames{i}]
            'BINARY'
            'O.K.'
            'DARK CURRENT'
            'YES'
            'DC FILE'
            DarkName
            'O.K.'
            'BINARY'
            'CAKE'
            'KEYBOARD'
            num2str(centerX, 7)
            num2str(centerY, 7)
            '1'
            '3.3627234E+03'
            '1.5874724E+03'
            '1'
            '3.3540134E+03'
            '1.4568239E+03'
            '1'
            '2.0649460E+03'
            '1.6223123E+03'
            '1'
            '3.3888530E+03'
            '1.3784344E+03'
            'INTEGRATE'
            'X-PIXEL SIZE'
            num2str(pixSize, 7)
            'Y-PIXEL SIZE'
            num2str(pixSize, 7)
            'DISTANCE'
            num2str(Dsam, 7)
            'WAVELENGTH'
            num2str(Lambda, 7)
            'TILT ROTATION'
            num2str(TiltPlane, 7)
            'ANGLE OF TILT'
            num2str(InPlane, 7)
            'O.K.'
            'START AZIMUTH'
            num2str(etaStart, 7)
            'END AZIMUTH'
            num2str(etaEnd, 7)
            'INNER RADIUS'
            num2str(rhoStart, 7)
            'OUTER RADIUS'
            num2str(rhoEnd, 7)
            'SCAN TYPE'
            'RADIAL'
            'AZIMUTH BINS'
            num2str(numBinsA)
            'RADIAL BINS'
            num2str(numBinsR)
            'CONSERVE INT.'
            'NO'
            'POLARISATION'
            'YES'
            'GEOMETRY COR.'
            'YES'
            'O.K.'
            'EXIT'
            'OUTPUT'
            'SPREAD SHEET'
            'YES'
            [SPRNames{i}]
            };
    else
        body{i} = {...
            'INPUT'
            [ImageNames{i}]
            'BINARY'
            'O.K.'
            'DARK CURRENT'
            'YES'
            'DC FILE'
            DarkName
            'O.K.'
            'BINARY'
            'CAKE'
            'INTEGRATE'
            'X-PIXEL SIZE'
            num2str(pixSize, 7)
            'Y-PIXEL SIZE'
            num2str(pixSize, 7)
            'X-BEAM CENTRE'
            num2str(centerX, 7)
            'Y-BEAM CENTRE'
            num2str(centerY, 7)
            'DISTANCE'
            num2str(Dsam, 7)
            'WAVELENGTH'
            num2str(Lambda, 7)
            'TILT ROTATION'
            num2str(TiltPlane, 7)
            'ANGLE OF TILT'
            num2str(InPlane, 7)
            'O.K.'
            'START AZIMUTH'
            num2str(etaStart, 7)
            'END AZIMUTH'
            num2str(etaEnd, 7)
            'INNER RADIUS'
            num2str(rhoStart, 7)
            'OUTER RADIUS'
            num2str(rhoEnd, 7)
            'SCAN TYPE'
            'RADIAL'
            'AZIMUTH BINS'
            num2str(numBinsA)
            'RADIAL BINS'
            num2str(numBinsR)
            'CONSERVE INT.'
            'NO'
            'POLARISATION'
            'YES'
            'GEOMETRY COR.'
            'YES'
            'O.K.'
            'EXIT'
            'OUTPUT'
            'SPREAD SHEET'
            'YES'
            [SPRNames{i}]
            };
    end
    % write body
    for j = 1:length(body{i})
        fprintf(fid, '%s\n', body{i}{j});
    end
end

trailer = {...
    'EXIT'
    'MACROS / LOG FILE'
    'CLOSE LOG FILE'
    'EXIT'
    'EXIT FIT2D'
    'YES'
    '%!*\ END OF IO MACRO FILE'
    };

for i = 1:length(trailer)
    fprintf(fid, '%s\n', trailer{i});
end

fclose(fid);