function q_pf = GeneratePFScattVectors(th, ome, eta, varargin)
% GeneratePFScattVectors - Generates pole figure scattering vectors for
% powder diffraction
%
%   USAGE:
%
%   qpf = GeneratePFScattVectors(th, ome, eta)
%
%   INPUT:
%
%   th
%       Bragg angle of the crystallographic plane of interest (scalar)
%       (deg)
%
%   ome
%       Rotation angles about YL (1 x n matrix where n is the number of
%       omega intervals) (deg)
%
%   eta
%       Azimuthal angles on the detector (1 x m matrix where m is the
%       number of eta intervals) (deg); refer to the following table for
%       angle convention.
%       eta = 0 => -X
%       eta = 90 => +Y
%       eta = 180 => +X
%       eta = 270 => -Y
%
%   Version (optional)
%       'default' is the version that existed before 2017-12-13
%       'aps' is the version that generates angles based on APS coordinate
%       system
%
%   OUTPUT:  none
%
%   qpf
%       a list of unit vectors on pole figure - 3 x (m x n)
%

% default options
optcell = {...
    'Version', 'default', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

ome = ome(:);
eta = eta(:);

n_ome   = length(ome);
n_eta   = length(eta);

switch lower(opts.Version)
    case 'default'
        % CONVERSION BETWEEN FIT2D +X AND APS +X
        eta = eta + 180;
        
        q0  = [1 0 0]';
        R_th    = [ ...
            cosd(th) 0 sind(th); ...
            0 1 0; ...
            -sind(th) 0 cosd(th); ...
            ];
        q0  = R_th*q0;

        q_eta   = zeros(3,n_eta);
        for iii = 1:1:n_eta
            R_eta   = [ ...
                cosd(eta(iii)) -sind(eta(iii)) 0; ...
                sind(eta(iii)) cosd(eta(iii)) 0; ...
                0 0 1; ...
                ]';
            q_eta(:,iii)  = R_eta*q0;
        end
        
        q_pf    = zeros(3,n_eta * n_ome);
        for iii = 1:1:n_ome
            R_ome   = [ ...
                cosd(ome(iii)) 0 sind(ome(iii)); ...
                0 1 0; ...
                -sind(ome(iii)) 0 cosd(ome(iii)); ...
                ];
            q_ome	= R_ome*q_eta;
            
            ini = 1 + n_eta * (iii - 1);
            fin = n_eta * iii;
            q_pf(:,ini:fin)  = q_ome;
        end
    case 'aps'
        q0  = [-1 0 0]';
        
        R_th    = [ ...
            cosd(th(1)) 0 sind(th(1)); ...
            0 1 0; ...
            -sind(th(1)) 0 cosd(th(1)); ...
            ];
        q1  = R_th*q0;
        
        q_eta   = zeros(3,n_eta);
        for iii = 1:1:n_eta
            Rz  = [ ...
                cosd(eta(iii)) -sind(eta(iii)) 0; ...
                sind(eta(iii)) cos(eta(iii)) 0; ...
                0 0 1]';
            q_eta(:,iii)    = Rz*q1;
        end
        
        q_pf    = zeros(3,n_eta * n_ome);
        for iii = 1:1:n_ome
            R_ome   = [ ...
                cosd(ome(iii)) 0 sind(ome(iii)); ...
                0 1 0; ...
                -sind(ome(iii)) 0 cosd(ome(iii)); ...
                ];
            q_ome	= R_ome*q_eta;
            
            ini = 1 + n_eta * (iii - 1);
            fin = n_eta * iii;
            q_pf(:,ini:fin)  = q_ome;
        end
    otherwise
        disp('convention not implemented')
end

% Data    = ones(size(q_pf,2),1);
% PlotSPF(q_pf', Data)