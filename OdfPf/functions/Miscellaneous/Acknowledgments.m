function a = Acknowledgments()
% Acknowledgments - Print acknowledgments.
%   
%   USAGE:
%
%   Acknowledgments
%   a = Acknowledgments
%
%   INPUT:  none
%
%   OUTPUT:
%
%   a is a cell array of strings,
%     containing the acknowledgments
%
fname = 'Acknowledgments.txt';
%
a = textread(fname, '%s', -1,   ... 
	     'delimiter', '\n', ...
	     'whitespace', '');
disp(strvcat(a))
