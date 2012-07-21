function header = ExportDXGroup(fname, gname, members, type)
% EXPORTDXGROUP - Export DX group.
%   
%   USAGE:
%
%   header = ExportDXGroup(fname, gname, members)
%
%   INPUT:
%
%   fname   is a string,
%           the full filename to write 
%   gname   is a string,
%           the name of the DX group object
%   members is a cell array of strings,
%           the list of group member names
%   type    is a string (optional)
%           specifies the group type, e.g. 'multigrid'
%
%   OUTPUT:
%
%   header is a cell array of strings,
%          the text of the DX group data file
%
%   NOTES:
%
%   * 'fname' is the full filename including suffix, not
%     the basename, as in 'ExportDX'
%   * the member filenames are assumed to be the membername
%     with .dx suffix
%
if nargin > 3
  gtype = type;
else
  gtype = 'group';
end
gclass = [' class ', gtype];
%
NL = sprintf('\n');
header = {['object ', quote(gname),  gclass]};
for i=1:length(members)
  member_i = members{i};
  memstring = ['   member ', quote(member_i), NL, ...
	       '      file ', member_i, '.dx'];
  header = {header{:}, memstring};
end
%
WriteDataFile(fname, header, []);
%
%-------------------- *
%
function qs = quote(s)
% quote - Enclose string in quotes.
%   
qs = ['"', s, '"'];
