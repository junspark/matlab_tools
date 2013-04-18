function [tth_range, d_range, Q_range] = vff_MaxCoverage(r_inner, r_outer, Dsam, lambda, varargin)
% vff_MaxCoverage - calculates the coverage of the very far field setup
%
%   USAGE:
%
%   vff_MaxCoverage(r_inner, r_outer, D, lambda)
%
%   INPUT:
%
%   r_inner
%       inner radius of the vff detector
%
%   r_outer
%       outer radius of the vff detector
%
%   Dsam
%       sample to detector distance
%
%   lambda
%       incident x-ray wavelength
%
%   OUTPUT:
%
%   tth_range
%       2-theta range
%
%   d_range
%       d-spacing range
%
%   Q_range
%       Q range
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'ScreenOutput'	prints out coverage results to screen

% default options
optcell = {...
    'ScreenOutput', 'off', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

tth_inner   = atand(r_inner/Dsam);
tth_outer   = atand(r_outer/Dsam);

d_inner = lambda/2/sind(tth_inner/2);
d_outer = lambda/2/sind(tth_outer/2);

Q_inner = 4*pi*sind(tth_inner/2)/lambda;
Q_outer = 4*pi*sind(tth_outer/2)/lambda;

tth_range   = [tth_inner; tth_outer];
d_range     = [d_inner; d_outer];
Q_range     = [Q_inner; Q_outer];

if strcmpi(opts.ScreenOutput, 'on')
    disp('************* COVERAGE INFO *************')
    disp(sprintf('  2theta_inner = %5.4f deg', tth_inner))
    disp(sprintf('  2theta_outer = %5.4f deg\n', tth_outer))
    disp(sprintf('  d_inner = %5.4f A @ %f Angstroms', d_inner, lambda))
    disp(sprintf('  d_outer = %5.4f A @ %f Angstroms\n', d_outer, lambda))
    disp(sprintf('  Q_inner = %5.4f A^-1 @ %f Angstroms', Q_inner, lambda))
    disp(sprintf('  Q_outer = %5.4f A^-1 @ %f Angstroms', Q_outer, lambda))
    disp('******************************************')
end
