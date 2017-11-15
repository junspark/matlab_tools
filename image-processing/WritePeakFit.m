function [] = WritePeakFit(pfname, pkfit, polimg, pklbl, tth, d)
% WritePkFit - writes out pkfit into a mat file
%
%   INPUT:
%
%   pfname
%           full file name of the output pkfit file
% 
%   pkfit
%       is the caked data generated using instr and cakeParms
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
save(pfname, 'pkfit', 'polimg', 'pklbl', 'tth', 'd')