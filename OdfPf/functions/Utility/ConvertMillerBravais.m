function outvec = ConvertMillerBravais(invec, varargin)
% CONVERTMILLERBRAVAIS - Convert array of Miller-Bravais vectors to/from
%   contracted notation: [hkil] <--> [hkl] where h + k + i = 0
%
% USAGE
%     outvec = ConvertMillerBravais(invec, options)
%
% INPUTS
%     1) invec is 3 x n or 4 x n, an array of Miller-Bravais vectors in
%          contracted or full notation, respectively
%     2) options is 1 x 1, a string containing either 'expand' or
%          'contract'.  In the case where contracted vectors are given and
%          'contract' is passed or full vectors are given and 'expand' is
%          passed, the output is set to the input
%
% OUTVEC
%     1) outvec is 3 x n or 4 x n, an array of Miller-Bravais vectors in
%          contracted or full notation, depending on the requested
%          functionality in options.
%   
% SEE ALSO
%     MillerBravaisToUnit
%
% NOTES
%     11/20/2003 *) There is a soft check of dimensionality compatability
%                     if full Miller-Bravais vectors are given, as well as
%                     between the input and the requested functionality.
n = size(invec, 1);
m = size(invec, 2);
%
if (isempty(varargin) & n == 3) | (strcmp(varargin{1}, 'expand') & n == 3)
    outvec = [invec(1, :);
              invec(2, :);
              -(invec(1, :) + invec(2, :));
              invec(3, :)];
elseif (isempty(varargin) & n == 4) | (strcmp(varargin{1}, 'contract') & n == 4)
    if sum(sum(invec(1:3, :), 1)) ~= 0
        error('Input array contains non-Miller-Bravais vectors')
    end
    outvec = [invec(1, :);
              invec(2, :);
              invec(4, :)];
elseif (strcmp(varargin{1}, 'expand') & n == 4)
    if sum(sum(invec(1:3, :), 1)) ~= 0
        error('Input array contains non-Miller-Bravais vectors')
    end
    outvec = invec;
elseif (strcmp(varargin{1}, 'contract') & n == 3)
    outvec = invec;
else
    error('Unkonwn/incompatable input args');
end