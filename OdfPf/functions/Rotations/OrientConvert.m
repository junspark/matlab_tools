function out = OrientConvert(in, inConv, outConv, inDeg, outDeg)
%  OrientConvert - Convert among orientation conventions.
%   
%   STATUS:  in development
%
%   USAGE:
%
%   out = OrientConvert(in, inConv, outConv)
%   out = OrientConvert(in, inConv, outConv, inDeg, outDeg)
%
%   INPUT:
%
%   in      is d x n 
%           input parameters (e.g. Euler angles)
%   inConv  is a string
%           input convention
%   outConv is a string
%           output convention
%   inDeg   is a string
%           either 'degrees' or 'radians'
%   outDeg  is a string
%           either 'degrees' or 'radians'
%
%   OUTPUT:
%
%   out is e x n 
%          output parameters
%   NOTES:
%
%   * Conventions are 'kocks', 'bunge', 'rmat', 'quat', 'rod'
%   * If any Euler angle conventions are specified, then the
%     degrees convention must also be specified
%
%

%  Convert input to rotation matrices.

switch inConv
 case 'kocks'
  rmat = RMatOfBunge(BungeOfKocks(in, inDeg), inDeg);
 case 'bunge'
  rmat = RMatOfBunge(in, inDeg);
 case 'rmat'
  rmat = in;
 case {'rod', 'rodrigues'}
  rmat = RMatOfQuat(QuatOfRod(in));
 case {'quat', 'quaternion'}
  rmat = RMatOfQuat(in);
 otherwise
  error('input convention not matched')
  
end

%  Convert rotation matrices to output.

switch outConv
 case 'kocks'
  out = KocksOfBunge(BungeOfRMat(rmat, outDeg), outDeg);
 case 'bunge'
  out = BungeOfRMat(rmat, outDeg);
 case 'rmat'
  out = rmat;
 case {'quat', 'quaternion'}
  out = QuatOfRMat(rmat);
 case {'rod', 'rodrigues'}
  out = RodOfQuat(QuatOfRMat(rmat));
 otherwise
  error('output convention not matched')
  
end
