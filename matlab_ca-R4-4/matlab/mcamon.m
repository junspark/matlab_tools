function varargout = mcamon(handle,varargin)
%MCAMON       - install or replace monitor on a PV
%   
% STS = MCAMON(HANDLE) installs monitor with default callback.
%   Default callback updates local copy of the channel data
%   every time the data changes on the server.
%   This cached data can be read at later time into MATLAB with MCACACHE.
%   HANDLE - integer handle to a channel previously opened with MCAOPEN
%   Returns 1 on success, 0 on failure
% 
% STS = MCAMON(HANDLE,CALLBACKSTRING) installs a monitor and specifies 
%   a callback string for each. A callback string must be a MATLAB command,
%   sequence of commands or a name of a script/function on the MATLAB path.
%   It is executed in the 'base' workspace (AFTER the default callback) on
%   the next poll of the queue by the MCAMONTIMER command.
%   Returns 1 on success, 0 on failure
%
% [HANDLES, CALLBACKSTRINGS]=MCAMON with no arguments returns information
%    on all currently installed monitors
%
% Note: Monitors can be installed with MCAMON and cleared with 
%    MCACLEARMON any number of times.  Use MCAMONTIMER to initialise
%    the MATLAB timer which polls and processes the outstanding MCA Monitor
%    callback queue.
% 
% Note:  Use of asynchronous features of EPICS (such as monitors) 
%    with MATLAB requires special care - read MCA notes.
%   
%    1.In CA client library (EPICS R3.13.4) asynchronous callbacks run one at a time
%    to completion. This means that MATLAB callback string installed with MCAMON
%    may not itself contain other MCA functions that call CA library such as MCAGET
%    For example MCAMON(H1, 'X=MCAGET(H2);') will not work.
%    MCAMON(H1, X='MCACACHE(H2);') is OK since MCACACHE does not use CA library. 
%
% See also MCAMONTIMER, MCACACHE, MCAGET, MCACLEARMON.

if nargin ==1
    varargout{1}=mca(100,handle);
elseif nargin==2 
    if ischar(varargin{1})
        varargout{1} = mca(100,handle,varargin{1});
    else
        error('Second argument must be a string');
    end
elseif nargin == 0
    if nargout == 2
        [varargout{1},varargout{2}]=mca(500);
    else
        varargout{1}=mca(500);    
    end
end