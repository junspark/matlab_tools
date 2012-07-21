function info = OdfPfVersion()
% OdfPfVersion - Print version information.
%   
%   USAGE:
%
%   OdfPfVersion;
%   info = OdfPfVersion;
%
%   INPUT:  none
%
%   OUTPUT:
%
%   info is a cell array of strings containing version information.
%
Y = '2007';
M = '01';
D = '08';
%
vname = sprintf('OdfPf-%s-%s-%s', Y, M, D);
vdate = sprintf('%s-%s-%s', Y, M, D);
%
svn = '04';
%
info = {...
    ['Version:  ', vname], ...
    ['   Date:  ', vdate], ...
    ['    Svn:  ', svn  ]  ...
       };
fprintf('%s\n', info{:});
