function [F, J] = sumfun1(pfree,pfix,npar,fnames,x,y,w)
% [F, J] = sumfun1(pfree,pfix,npar,fnames,x,y,w)
% Sums functions of common independent variable (multiple component LSQ fitting)
% INPUT
% pfree: FLOAT(1 < SUM(npar)), vector of free parameters to be varried in LSQ fiting 
% pfix: FLOAT(SUM(npar)), vector of fixed parameters separated by NaN
% npar: INT(nfuncs), vector containing the number of parameters of the functions in fnames
% fnames: CELL ARRAY, names of functions to be summed.
% x: FLOAT(npts), independent variable
% y: FLOAT(npts), dependent variable (when used for LSQ fitting)
% w: FLOAT(npts), weights for LSQ fitting
% OUTPUT
% F: FLOAT(npts), 5 input arguments: sum of function evaluation, 7 input arguments: (sum-y)/w
% J: FLOAT(npts, length(pfree)), partial derivatives
% COMMENTS
% Call from LSQNONLIN:
% [pfit,...]=LSQNONLIN('sumfun1',pfree,LB,UB,OPTIONS,pfix,npar,fnames,x,y,w)

if ~(nargin == 5)&&~(nargin == 7); error('5 or 7 input arguments must be supplied'); end
if isempty(pfix); pfix = ones(sum(npar),1)*NaN; end
if ~(length(pfix) == sum(npar)); error('length of pfix does not match number of parameters given by sum(npar)'); end
if ~(length(pfree) == sum(isnan(pfix))); error('length(pfree) does not match number of NaNs in pfix'); end;
if ~(length(npar) == length(fnames)); error('length(npar must be equal length(fnames))'); end
pfree = pfree(:); pfix = pfix(:); npar = npar(:); x = x(:);
nfun = length(npar);
p= pfix; p(isnan(pfix))=pfree;

pi = 1;
F = zeros(size(x));

if nargout == 1 % calculation of J not required
    for i = 1:nfun % calculate F
        F = feval(char(fnames(i)),p(pi:pi-1+npar(i)),x) + F;
        pi = pi + npar(i);
    end
    if nargin == 7 % if y & w are supplied, F= minimization function for LSQ fitting
        y = y(:); w = w(:);
        F = (F - y)./w;
    end
else % for 2 output arguments calculate F & J
    jn = 1; fi = zeros(size(x));
    for i = 1:nfun
        pind = pi:pi-1+npar(i);
        [fi, jp] = feval(char(fnames(i)),p(pind),x);
        F = F + fi;
        pi = pi + npar(i);
        J(:,jn:jn+sum(isnan(pfix(pind)))-1) = jp( :, isnan(pfix(pind))  );
        jn = jn + sum(isnan(pfix(pind)));
    end
    if nargin == 7 % if y & w are supplied, F= minimization function for LSQ fitting
        y = y(:); w = w(:);
        F = (F - y)./w;
        J = J./(w*ones(1,length(pfree)));
    end
end
   
   