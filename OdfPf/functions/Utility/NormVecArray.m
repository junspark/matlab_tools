function nrm = NormVecArray(mat, varargin)
% NORMVECARRAY - Generate the row or column 2-norm of an array of vectors
%
% USAGE
%      nrm = NormVecArray(mat, option)
%
% INPUTS
%     1) mat is n x m, an array of m n-vectors or n m-vectors
%     2) option is 1 x 1, a string.  By default, the column norm is taken.
%          To take the row norm, set option to 'rows'.
% OUTPUTS
%     1) nrm is eiter 1 x m or n x 1, a vector containing the column or row
%          norms depending on option. 
%
% SEE ALSO
%     norm, IPNorm

if length(varargin) == 1 & strcmp(varargin{1}, 'rows')
    dim = size(mat, 2);
    ddim = 2;
elseif isempty(varargin)
    dim = size(mat, 1);
    ddim = 1;
else
    error('Unknown input argument')
end

nrm = sqrt(dot(mat, mat, ddim));