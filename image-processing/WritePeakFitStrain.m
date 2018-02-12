function [] = WritePeakFitStrain(pfname_pkstrain, pfname_pkfit, ...
        strain)
% WritePeakFitStrain - writes out pkfit and associated strain into a mat file
%
%   INPUT:
%   pfname_pkstrain
%       full file path and name to the fit and strain information
%
%   pfname_pkfit
%       full file path and name to the fit information
% 
%   strain
%       calculated lattice strain
% 
%   OUTPUT:  
%           
%       none
save(pfname_pkstrain, 'pfname_pkfit', 'strain')