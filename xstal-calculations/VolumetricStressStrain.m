function tout = VolumetricStressStrain(tin)
% VolumetricStressStrain - calculates the volumetric stress or strain
%
%   USAGE:
%
%   tout = VolumetricStressStrain(tin)
%
%   INPUT:
%
%   tin
%       modified Voight-Mandel stress (strain) vector in [11, 22, 33, 12,
%       13, 23] order. Shears are expected to have factor of sqrt(2)
%       multiplied to them already.
%
%   tout
%       volumetric / hydrostatic stress (strain)

tin = tin(:);
if length(tin) == 6
    h   = sum(tin(1:3))/3;
    tout    = h;
else
    disp('input stress vector size incorrect');
end