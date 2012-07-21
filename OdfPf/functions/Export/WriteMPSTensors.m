function h = WriteMPSTensors(fname, t)
% WriteMPSTensors - Write tensors in MPS format.
%   
%   STATUS:  in development
%
%   USAGE:
%
%   h = WriteMPSTensors(t)
%
%   INPUT:
%
%   fname is a string
%         the name of the output file
%   t     is m x n
%         an array of n m-vectors
%
%   OUTPUT:
%
%   h is a cell array of strings
%        the MPS header of the data file
%
%   NOTES:
%
%   * For general tensor, expects t to be 9 x n (not 3 x 3 x n)
%
switch size(t, 1)
 case  1
  ttype = 'spherical';
 case  3
  ttype = 'skew'; 
 case  5
  ttype = 'deviatoric';
 case  6
  ttype = 'symmetric';
 case  9
  ttype = 'general';
 otherwise
  error('data first dimension does not correspond to any tensor format')
end
%
ntens = size(t, 2);
%
h = {...
    sprintf('%%number-of-tensors %d', ntens), 
    sprintf('%%type-of-tensor    %s', ttype),
    '%end-header'};
%
WriteDataFile(fname, h, t');
