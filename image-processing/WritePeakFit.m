function [] = WritePeakFit(pfname, ...
        amp_fit, fwhm_fit, mix_fit, tth_fit, bkg_fit, rsn_fit, ef_fit, rwp_fit, ...
        polimg, pklbl, tth, d)
% WritePkFit - writes out pkfit into a mat file
%
%   INPUT:
%
%   pfname
%           full file name of the output pkfit file
% 
%   amp, fwhm, mix, tth, bkg, rsn, ef, rwp
%       are the peak fit information
% 
%   polimg
%       is the caked data generated using instr and cakeParms
%
%   pklbl
%       peak labels
%
%   tth
%      2theta at reference state
%
%   d
%      d-spacings at the reference state
%
%   OUTPUT:  
%           
%           none

pkfit.amp_fit   = amp_fit;
pkfit.fwhm_fit  = fwhm_fit;
pkfit.mix_fit   = mix_fit;
pkfit.tth_fit   = tth_fit;
pkfit.bkg_fit   = bkg_fit;
pkfit.rsn_fit   = rsn_fit;
pkfit.ef_fit    = ef_fit;
pkfit.rwp_fit   = rwp_fit;
            
save(pfname, 'pkfit', 'polimg', 'pklbl', 'tth', 'd')