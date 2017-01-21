function tout = DeviatoricStressStrain(tin)
% DeviatoricStressStrain - calculates the deviatoric stress or strain
%
%   USAGE:
%
%   tout = DeviatoricStressStrain(tin)
%
%   INPUT:
%
%   tin
%       modified Voight-Mandel stress (strain) vector in [11, 22, 33, 12,
%       13, 23] order. Shears are expected to have factor of sqrt(2)
%       multiplied to them already.
%
%   tout
%       deviatoric stress (strain). the output is [11, 22, 33, 12, 13, 23]
%       order. Shears are expected to have factor of sqrt(2) multiplied to
%       them

tin = tin(:);
if length(tin) == 6
    h           = VolumetricStressStrain(tin);
    tout        = tin;
    tout(1:3)   = tin(1:3) - h;
else
    disp('input stress vector size incorrect');
end