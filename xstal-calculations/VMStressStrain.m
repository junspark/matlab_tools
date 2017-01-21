function vm = VMStressStrain(tin)
% VMStressStrain - calculates the von Mises effective stress or strain
%
%   USAGE:
%
%   vm = VMStressStrain(tin)
%
%   INPUT:
%
%   tin
%       modified Voight-Mandel stress (strain) vector in [11, 22, 33, 12,
%       13, 23] order. Shears are expected to have factor of sqrt(2)
%       multiplied to them already.
%
%   tout
%       von Mises effective stress (strain).

tin = tin(:);
if length(tin) == 6
    vm  = (tin(1) - tin(2))^2 + (tin(2) - tin(3))^2 + (tin(3) - tin(1))^2 + 3*(tin(4)^2 + tin(5)^2 + tin(6)^2);
    vm  = sqrt(vm/2);
else
    disp('input vectror size incorrect');
end