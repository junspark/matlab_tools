function Material = GSAS2_MakeMaterialFromPkslst(pfname_pkslst, pfname_imctrl, varargin)
% GSAS2_MakeMaterialFromPkslst - Makes material from GSAS2 PWDR pkslst
%
%   INPUT:
%
%   pfname_pkslst
%       full pathname to the GSAS2 pkslst file.
%   
%   pname_chi
%       pathname to where the chi files exist.
%
%   froot_chi
%       file name root of rht chi files.
%
%   OUTPUT:
%
%   Material
%       material compatible with FitPeaksPerPolImage 

% default options
optcell = {...
    'd_spacing_range', 0.03, ...
    };

% update option
opts            = OptArgs(optcell, varargin);
d_spacing_range = opts.d_spacing_range;
d_spacing_range = d_spacing_range(:);

% clear all
% close all
% clc
% 
% pkslst  = GSAS2_Read_pkslst('D:\w\mpe_apr21_analysis\park_sae_phi0_2\mpe_apr21_park_sae_ge3_80keV_2278mm_phi0.pkslst');
% imctrl  = GSAS2_Read_imctrl('D:\w\mpe_apr21_analysis\park_sae_phi0_2\mpe_apr21_park_sae_ge3_80keV_2278mm_phi0.imctrl');
% % d_spacing_range = 0.03;
% d_spacing_range = [0.02 0.02 0.015 0.015 0.015 0.015]';

%%% LOAD GSAS2 GENERATED FILES
pkslst  = GSAS2_Read_pkslst(pfname_pkslst);
imctrl  = GSAS2_Read_imctrl(pfname_imctrl);

%%% MATERIAL CALCULATION
numpks      = size(pkslst,1);
numbounds   = length(d_spacing_range);
tth         = pkslst(:,1);
th          = tth./2;
d           = (imctrl.wavelength/2)./sind(th);

if numbounds == 1
    d_spacing_UB    = (1 + d_spacing_range).*d;
    d_spacing_LB    = (1 - d_spacing_range).*d;
elseif numbounds == numpks
    d_spacing_UB    = (1 + d_spacing_range).*d;
    d_spacing_LB    = (1 - d_spacing_range).*d;
else
    disp('something odd with numpks')
    Material    = [];
    return
end
tth_UB  = 2.*asind(imctrl.wavelength/2)./d_spacing_LB;
tth_LB  = 2.*asind(imctrl.wavelength/2)./d_spacing_UB;

%%% MAKE MATERIALS VARIABLE
Material.num        = 1;
Material.pkbck      = 2;
Material.pkfunc     = 4;

% Material.lattparms  = 3.60;
% Material.structure  = 'fcc';
Material.numbounds  = numbounds;
Material.numpk      = numpks;
Material.tth        = tth;
Material.tth_UB     = tth_UB;
Material.tth_LB     = tth_LB;
Material.d_spacing  = d;

for iii = 1:1:numpks
    Material.pkidx{iii,1}   = iii;
end

Material.pkrange    = [
    tth_LB'; ...
    tth_UB'; ...
    ];

