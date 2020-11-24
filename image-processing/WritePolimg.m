function [] = WritePolimg(pfname, polimg, instr, omega, chi, cakeprms)
% WritePolimg - writes out polimg into a mat file
%
%   INPUT:
%
%   pfname
%           full file name of the output polimg file
% 
%   polImg
%           is the caked data generated using instr and cakeParms
% 
%   instr
%           is the instrument parameters in a structure array
%
%   omega
%           is the rotation about the Y axis in APS
%
%   chi
%           is the rotation about the Z axis in APS (rarely used)
%
%   cakeParms
%           is the structure array that defines the caking parameters
%
%   OUTPUT:  
%           
%           none

instr.omega = omega;
instr.chi   = chi;

save(pfname, 'polimg', 'instr', 'cakeprms')
