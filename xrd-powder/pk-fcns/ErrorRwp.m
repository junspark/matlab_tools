function Rwp = ErrorRwp(yData, yCalc)
% ErrorRwp - calculates weighted residual of the pattern.  for more details
% look at Rietvald Method (Young).
%
%   USAGE:
%
%   Rwp = ErrorRwp(yData, yCalc)
%
%   INPUT:
%   yData
%       experimental intensity data
%
%   yCalc
%       calculated intensity
%
%   OUTPUT:
%
%   Rwp
%       weighted residual of the pattern

i   = yData == 0;
yData(i)    = eps;

w   = 1./abs(yData);
Rwp = sum(w.*((yData - yCalc).^2))/sum(w.*(yData.^2));
Rwp = Rwp.^.5;