function q_pf = GeneratePFScattVectors(th, ome, eta)
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
%       Rotation angles about YS (1 x n matrix where n is the number of
%       omega intervals) (deg)
%
%   eta
%       Azimuthal angles on the detector (1 x m matrix where m is the
%       number of eta intervals) (deg)
%
%   OUTPUT:  none
%
%   qpf
%       a list of unit vectors on pole figure - 3 x (m x n)
%

% CONVERSION BETWEEN FIT2D +X AND APS +X
eta = eta + 180;

q0  = [1 0 0]';
R_th    = [ ...
    cosd(th) 0 sind(th); ...
    0 1 0; ...
    -sind(th) 0 cosd(th); ...
    ];
q0  = R_th*q0;

n_ome   = length(ome);
n_eta   = length(eta);

q_eta   = zeros(3,n_eta);
for i = 1:1:n_eta
    R_eta   = [ ...
        cosd(eta(i)) -sind(eta(i)) 0; ...
        sind(eta(i)) cosd(eta(i)) 0; ...
        0 0 1; ...
        ]';
    q_eta(:,i)  = R_eta*q0;
end

q_pf    = zeros(3,n_eta * n_ome);
for i = 1:1:n_ome
    R_ome   = [ ...
        cosd(ome(i)) 0 sind(ome(i)); ...
        0 1 0; ...
        -sind(ome(i)) 0 cosd(ome(i)); ...
        ];
    q_ome	= R_ome*q_eta;
    
    ini = 1 + n_eta * (i - 1);
    fin = n_eta * i;
    q_pf(:,ini:fin)  = q_ome;
end

% Data    = ones(size(q_pf,2),1);
% PlotSPF(q_pf', Data)