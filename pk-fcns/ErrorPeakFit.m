function variance = ErrorPeakFit(pr, resd, jacobian)
%   ErrorPeakFit generates various goodness of peak fit error metrics
%
%   USAGE:
%   output = ErrorPeakFit()
%
%   INPUT:
%
%   pr
%   peak fit parameters
%
%   resd
%   residual from the peak fit
%
%   jacobian
%   jacobian of the peak fit
%
%

[conf, variance]    = confint(pr, resd, jacobian); % variance of fitted params
variance            = full(variance);

% Afit{ppp}(mmm,1)                = p(1);
% Afit_errorbar{ppp}(mmm,1)       = sqrt(var(1));
% Afit_conf95_range{ppp}(mmm,1)   = conf(1,2) - conf(1,1);
% 
% gfit{ppp}(mmm,1)                = p(2);
% gfit_errorbar{ppp}(mmm,1)       = sqrt(var(2));
% gfit_conf95_range{ppp}(mmm,1)   = conf(2,2) - conf(2,1);
% 
% nfit{ppp}(mmm,1)                = p(3);
% nfit_errorbar{ppp}(mmm,1)       = sqrt(var(3));
% nfit_conf95_range{ppp}(mmm,1)   = conf(3,2) - conf(3,1);
% 
% Efit{ppp}(mmm,1)                = p(4);
% Efit_errorbar{ppp}(mmm,1)       = sqrt(var(4));
% Efit_conf95_range{ppp}(mmm,1)   = conf(4,2) - conf(4,1);
% 
% bkg{ppp}{mmm,1}                 = p(5:end);
% bkg_errorbar{ppp}{mmm,1}        = sqrt(var(5:end));
% bkg_conf95_range{ppp}{mmm,1}    = conf(5:end,2) - conf(5:end,1);
% 
% Rwp{ppp}(mmm,1) = ErrorRwp(ydata, yfit);
% Re{ppp}(mmm,1)  = ErrorRe(ydata, yfit);
% Rp{ppp}(mmm,1)  = ErrorRp(ydata, yfit);
% Rn{ppp}(mmm,1)  = rn;