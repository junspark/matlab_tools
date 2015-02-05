clear all
close all
clc

%%%%%%%%%%%
CuNb_500nm_parallel = load('W:\__eval\park_nov14\ge\CuNb_500nm_09019.ge3.frame1.cor.polimg.mat');
Data    = cell(1, CuNb_500nm_parallel.XRDIMAGE.CakePrms.bins(1));
for ii=1:1:CuNb_500nm_parallel.XRDIMAGE.CakePrms.bins(1)
    Data{ii}    = [CuNb_500nm_parallel.XRDIMAGE.Instr.pixelsize*CuNb_500nm_parallel.polimg.radius(ii,:)' CuNb_500nm_parallel.polimg.intensity(ii,:)'];
end

%%% DEPENDS ON WHICH MODEL
mapped_tth  = GeometricModelXRD2a(...
    CuNb_500nm_parallel.XRDIMAGE.Instr.centers./1000, ...
    CuNb_500nm_parallel.XRDIMAGE.Instr.distance, ...
    CuNb_500nm_parallel.XRDIMAGE.Instr.gammaY, CuNb_500nm_parallel.XRDIMAGE.Instr.gammaX, ...
    Pixel2mm(CuNb_500nm_parallel.polimg.radius', CuNb_500nm_parallel.XRDIMAGE.Instr.pixelsize), CuNb_500nm_parallel.polimg.azimuth, CuNb_500nm_parallel.XRDIMAGE.Instr.detpars)';

CuNb_500nm_parallel.polimg.mapped_tth_for_intensity = mapped_tth;

[tth_grid, intensity_in_tth_grid]   = MapIntensityToTThGrid(CuNb_500nm_parallel.XRDIMAGE, CuNb_500nm_parallel.polimg);
CuNb_500nm_parallel.polimg.tth_grid                 = tth_grid;
CuNb_500nm_parallel.polimg.intensity_in_tth_grid    = intensity_in_tth_grid;

%%%%%%%%%%%
CuNb_500nm_norm     = load('W:\__eval\park_nov14\ge\CuNb_500nm_09019.ge3.frame91.cor.polimg.mat');
Data    = cell(1, CuNb_500nm_norm.XRDIMAGE.CakePrms.bins(1));
for ii=1:1:CuNb_500nm_norm.XRDIMAGE.CakePrms.bins(1)
    Data{ii}    = [CuNb_500nm_norm.XRDIMAGE.Instr.pixelsize*CuNb_500nm_norm.polimg.radius(ii,:)' CuNb_500nm_norm.polimg.intensity(ii,:)'];
end

%%% DEPENDS ON WHICH MODEL
mapped_tth  = GeometricModelXRD2a(...
    CuNb_500nm_norm.XRDIMAGE.Instr.centers./1000, ...
    CuNb_500nm_norm.XRDIMAGE.Instr.distance, ...
    CuNb_500nm_norm.XRDIMAGE.Instr.gammaY, CuNb_500nm_norm.XRDIMAGE.Instr.gammaX, ...
    Pixel2mm(CuNb_500nm_norm.polimg.radius', CuNb_500nm_norm.XRDIMAGE.Instr.pixelsize), CuNb_500nm_norm.polimg.azimuth, CuNb_500nm_norm.XRDIMAGE.Instr.detpars)';

CuNb_500nm_norm.polimg.mapped_tth_for_intensity = mapped_tth;

[tth_grid, intensity_in_tth_grid]   = MapIntensityToTThGrid(CuNb_500nm_norm.XRDIMAGE, CuNb_500nm_norm.polimg);
CuNb_500nm_norm.polimg.tth_grid                 = tth_grid;
CuNb_500nm_norm.polimg.intensity_in_tth_grid    = intensity_in_tth_grid;

%%%%%%%%%%%
CuNb_30nm_parallel  = load('W:\__eval\park_nov14\ge\CuNb_30nm_09041.ge3.frame91.cor.polimg.mat');
Data    = cell(1, CuNb_30nm_parallel.XRDIMAGE.CakePrms.bins(1));
for ii=1:1:CuNb_30nm_parallel.XRDIMAGE.CakePrms.bins(1)
    Data{ii}    = [CuNb_30nm_parallel.XRDIMAGE.Instr.pixelsize*CuNb_30nm_parallel.polimg.radius(ii,:)' CuNb_30nm_parallel.polimg.intensity(ii,:)'];
end

%%% DEPENDS ON WHICH MODEL
mapped_tth  = GeometricModelXRD2a(...
    CuNb_30nm_parallel.XRDIMAGE.Instr.centers./1000, ...
    CuNb_30nm_parallel.XRDIMAGE.Instr.distance, ...
    CuNb_30nm_parallel.XRDIMAGE.Instr.gammaY, CuNb_30nm_parallel.XRDIMAGE.Instr.gammaX, ...
    Pixel2mm(CuNb_30nm_parallel.polimg.radius', CuNb_30nm_parallel.XRDIMAGE.Instr.pixelsize), CuNb_30nm_parallel.polimg.azimuth, CuNb_30nm_parallel.XRDIMAGE.Instr.detpars)';

CuNb_30nm_parallel.polimg.mapped_tth_for_intensity = mapped_tth;

[tth_grid, intensity_in_tth_grid]   = MapIntensityToTThGrid(CuNb_30nm_parallel.XRDIMAGE, CuNb_30nm_parallel.polimg);
CuNb_30nm_parallel.polimg.tth_grid                 = tth_grid;
CuNb_30nm_parallel.polimg.intensity_in_tth_grid    = intensity_in_tth_grid;

%%%%%%%%%%%
CuNb_30nm_norm      = load('W:\__eval\park_nov14\ge\CuNb_30nm_09041.ge3.frame1.cor.polimg.mat');
Data    = cell(1, CuNb_30nm_norm.XRDIMAGE.CakePrms.bins(1));
for ii=1:1:CuNb_30nm_norm.XRDIMAGE.CakePrms.bins(1)
    Data{ii}    = [CuNb_30nm_norm.XRDIMAGE.Instr.pixelsize*CuNb_30nm_norm.polimg.radius(ii,:)' CuNb_30nm_norm.polimg.intensity(ii,:)'];
end

%%% DEPENDS ON WHICH MODEL
mapped_tth  = GeometricModelXRD2a(...
    CuNb_30nm_norm.XRDIMAGE.Instr.centers./1000, ...
    CuNb_30nm_norm.XRDIMAGE.Instr.distance, ...
    CuNb_30nm_norm.XRDIMAGE.Instr.gammaY, CuNb_30nm_norm.XRDIMAGE.Instr.gammaX, ...
    Pixel2mm(CuNb_30nm_norm.polimg.radius', CuNb_30nm_norm.XRDIMAGE.Instr.pixelsize), CuNb_30nm_norm.polimg.azimuth, CuNb_30nm_norm.XRDIMAGE.Instr.detpars)';

CuNb_30nm_norm.polimg.mapped_tth_for_intensity = mapped_tth;

[tth_grid, intensity_in_tth_grid]   = MapIntensityToTThGrid(CuNb_30nm_norm.XRDIMAGE, CuNb_30nm_norm.polimg);
CuNb_30nm_norm.polimg.tth_grid                 = tth_grid;
CuNb_30nm_norm.polimg.intensity_in_tth_grid    = intensity_in_tth_grid;
