clear all
close all
clc


%%% SX ELASTICITY TENSOR 
%%% SHEARS e12 NOT g12 THUS 2 MULTIPLIER
%%% STRESS AND STRAIN ARE VECTORIZED WITH SQRT(2) ON SHEAR TERMS

%%% HCP - Ti811
C11 = 141;
C33 = 163;
C12 = 76.9;
C13 = 57.9;
C44 = 48.7;
C66 = (C11 - C12)/2;
C   = [C11 C33 C12 C13 2*C44 C66];
S   = C2S(C, 'Structure', 'hexagonal');
S	= [ ...
    S(1) S(3) S(4) 0    0    0;
    S(3) S(1) S(4) 0    0    0;
    S(4) S(4) S(2) 0    0    0;
    0    0    0    S(6) 0    0;
    0    0    0    0    S(5) 0;
    0    0    0    0    0    S(5);
    ];

%%% HCP MAGNESIUM
% S   = [22.15 19.75 -7.7 -4.93 30.15 29.85]./1000;
% S	= [ ...
%     S(1) S(3) S(4) 0    0    0;
%     S(3) S(1) S(4) 0    0    0;
%     S(4) S(4) S(2) 0    0    0;
%     0    0    0    S(6) 0    0;
%     0    0    0    0    S(5) 0;
%     0    0    0    0    0    S(5);
%     ];


%%% CUBIC - Cu
% C   = [166 120 2*76];
% S   = C2S(C);
% S   = [...
%     S(1) S(2) S(2) 0    0    0;
%     S(2) S(1) S(2) 0    0    0;
%     S(2) S(2) S(1) 0    0    0;
%     0    0    0    S(3) 0    0;
%     0    0    0    0    S(3) 0;
%     0    0    0    0    0    S(3);
%     ];

%%% CUBIC - Al
% C   = [107 60.9 2*28.3];
% S   = C2S(C);
% S   = [...
%     S(1) S(2) S(2) 0    0    0;
%     S(2) S(1) S(2) 0    0    0;
%     S(2) S(2) S(1) 0    0    0;
%     0    0    0    S(3) 0    0;
%     0    0    0    0    S(3) 0;
%     0    0    0    0    0    S(3);
%     ];

%%% ISOTROPIC
% C   = [134.6 57.7 76.9];
% S   = C2S(C);
% S   = [...
%     S(1) S(2) S(2) 0    0    0;
%     S(2) S(1) S(2) 0    0    0;
%     S(2) S(2) S(1) 0    0    0;
%     0    0    0    S(3) 0    0;
%     0    0    0    0    S(3) 0;
%     0    0    0    0    0    S(3);
%     ];

%%%%%%%%%%%%%%%%%%%%%%
%%% FULL SPHERE
% n       = 40;
% [x,y,z] = sphere(n);
% x       = x(:);
% y       = y(:);
% z       = z(:);
% v       = [x y z];
% n2      = (n+1)^2;

%%% ONLY THE FUNDAMENTAL TRIANGLE
%%% HCP
ph1 = 0:1:90;
ph2 = 0:1:30;

%%% CUBIC
% ph1 = 45:1:90;
% ph2 = 0:1:35.26;

x   = cosd(ph1');
y   = zeros(length(ph1),1);
z   = sind(ph1');
v0  = [x y z];

v   = [];
for i = 1:1:length(ph2)
    Rz  = [...
        cosd(ph2(i)) -sind(ph2(i)) 0; ...
        sind(ph2(i)) cosd(ph2(i))  0; ...
        0            0             1; ...
        ];
    v   = [v Rz*v0'];
end
v   = v';
n2  = size(v,1);

%%%%%%%%%%%%%%%%%%%%%%
E_d     = zeros(n2,1);      % SAME UNITS AS INPUT PROPERTIES
E       = zeros(n2,1);      % SAME UNITS AS INPUT PROPERTIES
nuxy    = zeros(n2,1);
nuxz    = zeros(n2,1);
for i = 1:1:n2
    a   = v(i,:)';
    
    if a(3) == 1 | a(3) == -1
        PHI = -90;
        TH  = 0;
    elseif a(3) == 0
        PHI = 0;
        TH  = acosd(dot(a,[1 0 0]'));
    else
        axy = [a(1:2); 0];
        axy = axy./norm(axy);
        
        PHI = -acosd(dot(axy,a));
        TH  = acosd(dot(axy,[1 0 0]'));
    end
    
%     pause
%     TH  = -0;
%     PHI = 20.6986;
    
    Ry  = [ ...
        cosd(PHI)  0 sind(PHI); ...
        0          1 0        ; ...
        -sind(PHI) 0 cosd(PHI); ...
        ];
    Rz  = [ ...
        cosd(TH) -sind(TH) 0; ...
        sind(TH) cosd(TH)  0; ...
        0        0         1; ...
        ];
    
    R   = Rz*Ry;        % Rs = c
    T   = VectorizedCOBMatrix(R);
    
    S_d = T'*S*T;
    
    E_d(i)  = 1./S_d(1,1);
       
    nuxy(i) = -S_d(1,2)/S_d(1,1);
    nuxz(i) = -S_d(1,3)/S_d(1,1);
    
    %%% USE HOSFORD EQN TO VERIFY
    %%% CUBIC
    E(i)    = 1/(S(1,1) - 2*(S(1,1) - S(1,2) - S(4,4))*(a(1)*a(1)*a(2)*a(2) + a(2)*a(2)*a(3)*a(3) + a(3)*a(3)*a(1)*a(1)));
    
    %%% HCP
    g       = dot(a,[0 0 1]');
    E(i)    = ...
        1/(S(1,1) * (1-g^2)^2 + ...
        S(3,3) * g^4 + ...
        (g^2) * (1 - g^2) * (2*S(1,3) + 2*S(6,6)));
end
figure(1)
PlotSPF(v, E_d, ...
    'Title', 'E_{d} (MPa)', 'ViewAngle', [0 90], ...
    'ShowSurf', 'off')
figure(2)
PlotSPF(v, E, ...
    'Title', 'E_{d} (MPa)', 'ViewAngle', [0 90], ...
    'ShowSurf', 'off')
figure(3)
PlotSPF(v, nuxy, ...
    'Title', '\nu_{xy}', 'ViewAngle', [0 90], ...
    'ShowSurf', 'off')
figure(4)
PlotSPF(v, nuxz, ...
    'Title', '\nu_{xz}', 'ViewAngle', [0 90], ...
    'ShowSurf', 'off')