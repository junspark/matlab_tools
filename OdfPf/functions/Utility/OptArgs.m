function opts = OptArgs(optkeys, optargs)
% OptArgs - Set options from cell array.
%   
%   USAGE:
%
%   opts = OptArgs(optkeys, optargs)
%
%   INPUT:
%
%   optkeys is a cell array, (of length 2*n)
%           consisting of data in the form
%           {string1, object1, string2, object2, ...};
%           the strings are the parameter names and objects
%           are the default values associated with each paramter
%   optargs is a cell array, (of length 2*m)
%           it is of the same form as `optkeys'; its values
%           override the defaults 
%   
%   OUTPUT:
%
%   opts is a structure,
%           the fields are parameter names from `optkeys', and
%           the values are either the default value or the 
%           overriding value from `optargs'
%
%   NOTES:  
%
%   *  Values in `optkeys' which are cell arrays need to be
%      contained inside another cell array; see `struct' for
%      details, whereas values in `optargs' which are cell arrays 
%      should not be placed in another cell array.  The difference
%      arises because the resulting structure is generated by
%      a call to the matlab `struct' command with argument
%      `optkeys', but the structure fields are updated one by one
%      using the values in `optargs'.
%
%   *  The general usage of this function would be to set optkeys
%      in the user function and pass varargin as optargs to update
%      the default values.
%
opts    = struct(optkeys{:});
optargs = reshape(optargs, [2 length(optargs)/2]);
%
for k=1:size(optargs, 2)
  %
  key = optargs{1, k};
  val = optargs{2, k};
  %
  if (isfield(opts, key))
    opts.(key) = val;
  else
    warning('OptArgs:NoKey', ['no such option:  ', key])
  end
end