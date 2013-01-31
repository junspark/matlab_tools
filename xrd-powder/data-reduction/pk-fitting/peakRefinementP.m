function funcOut=peakRefinementP(x,xdata)

% A=x(1);
% mix=x(2);
% G=x(3);
% x0=x(4);
% bkG1=x(5);
% bkG2=x(6);
bkGrndI=xdata.*x(5)+x(6);
% spectra=bkGrndI;
% 
% A  = p(1);
% n  = p(2);
% G  = p(3);
% x0 = p(4);

% the TCH pseudo-Voight:
%
% f = A*((1 - n)*fg + n*fl)

% p=[A mix G
% pseudoVoight1([x(1) x(2) x(3) x(4)], xdata)
bkGrndI=xdata.*x(5)+x(6);
func = pseudoVoight1([x(1) x(2) x(3) x(4)], xdata);
spectra=bkGrndI+func';


funcOut=spectra;
