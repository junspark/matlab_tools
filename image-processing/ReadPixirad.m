function csq = ReadPixirad(pfname, varargin)
% ReadPixirad - read Pixirad hexagonal grid file
%
%   INPUT:
%
%   pfname
%       name of the Pixirad image file
%
%   nx
%       number of horizontal nodes (default = 512)
%
%   ny
%       number of vertical nodes (default = 476)
%
%   OUTPUT:
%
%   csq
%       Pixirad image file information mapped to an image with square
%       pixels
%

% default options
optcell = {...
    'nx', 512, ...
    'ny', 476, ...
    'nxsq', 512, ...
    'nysq', 476, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

% read in image
c           = imread(pfname);
[nx, ny]    = size(c);

if nx ~= opts.nx
    disp(sprintf('user input or default : %d', opts.nx))
    disp(sprintf('image size in x       : %d', nx))
    error('nx does not match')
elseif ny ~= opts.ny
    disp(sprintf('user input or default : %d', opts.ny))
    disp(sprintf('image size in y       : %d', ny))
    error('ny does not match')
else
    disp(sprintf('%d x %d image', nx, ny))
end
c 	= double(c(:));

pfname_xymap    = ['pixirad.map.nx.', num2str(nx), '.y.', num2str(ny), '.mat'];
if exist(pfname_xymap, 'file')
    load(pfname_xymap)
else
    disp(sprintf('pixirad map %s does not exist!', pfname_xymap))
    disp(sprintf('creating %s.', pfname_xymap))
    x1  = 1:1:nx;
    x2  = 0.5:1:(nx-0.5);
    
    x   = [];
    y   = [];
    
    ct  = 1;
    for i = 1:1:ny
        if mod(i,2) == 1
            x   = [x; x1'];
            y   = [y; sind(60)*ct*ones(length(x1),1)];
        else
            x   = [x; x2'];
            y   = [y; sind(60)*ct*ones(length(x2),1)];
        end
        ct  = ct + 1;
    end
    save(pfname_xymap, 'x', 'y')
end

%%% GENERATE SQUARE GRID
xsq = linspace(min(x), max(x), opts.nxsq);
ysq = linspace(min(y), max(y), opts.nysq);
[xsq, ysq]  = meshgrid(xsq, ysq);

%%% CREATE INTERPOLATION
F   = TriScatteredInterp(x, y, c);

%%% MAP HEX GRID DATA TO SQUARE GRID
csq = F(xsq, ysq);