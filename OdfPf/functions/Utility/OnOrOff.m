function toggle = OnOrOff(string)
% OnOrOff - Convert on|off string to 1|0 integer.
%
%   USAGE:
%
%   toggle = OnOrOff(string)
%
%   INPUT:
%
%   string is either 'on' or 'off' (case is ignored)
%
%   OUTPUT:
%
%   toggle is a scalar, 
%          0   if string matches 'off'
%          1   if string matches 'on'
%  
string = lower(string); % to lower case
%
if     (strcmp(string, 'on'))
  toggle = 1;
elseif (strcmp(string, 'off'))
  toggle = 0;
else
  error('string is not either ''on'' or ''off''')
end
