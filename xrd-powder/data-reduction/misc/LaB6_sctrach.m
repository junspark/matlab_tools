%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% MATERIAL PARAMETERS - LaB6
% XRDIMAGE.Material.num       = 1;
% XRDIMAGE.Material.lattparms = 4.1569162;        % LaB6
% XRDIMAGE.Material.structure = 'simplecubic';
% XRDIMAGE.Material.numpk     = 8;
% XRDIMAGE.Material.pkrange    = [...
%     2.7  3.8 4.7 6.1 6.7 8.2 8.7 9.1; ...
%     2.95 4.1 5.0 6.7 7.0 8.5 9.0 9.4; ...
%     ];
% XRDIMAGE.Material.pkidx     = {...
%     [1] [2] [3] [5] [6] [8] [9] [10]
%     };
% XRDIMAGE.Material.pkbck     = 2;
% XRDIMAGE.Material.pkfunc    = 4;
% XRDIMAGE.Material.hkls      = load([XRDIMAGE.Material.structure, '.hkls']);
% 
% %%% CALCULATE THEORETICAL TTH
% [d, th] = PlaneSpacings(XRDIMAGE.Material.lattparms, ...
%     'cubic', XRDIMAGE.Material.hkls', ...
%     XRDIMAGE.Instr.wavelength);
% tth     = 2*th;
% 
% XRDIMAGE.Material.tth       = tth;
% XRDIMAGE.Material.d_spacing = d;