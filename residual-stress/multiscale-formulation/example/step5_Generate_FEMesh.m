% GENERATE MESH
% Use MeshandBCs/DiskMesh.m to generate mesh for 3D solver
% Need to edit input.m to generate appropriate mesh
% Output is meshdata.mat
% Move meshdata.mat to Residual_Stress_3D_MLSEL_v1/ for step5

clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

PFNAME_GRID = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_GRID);
load(PFNAME_GRID)

%%% FIND NODAL POINTS
x	= linspace(0, 1.8, 3);
y	= linspace(0, 1.8, 3);
z	= linspace(0, 3, 6);

nox = size(x,2);
noy = size(y,2);
noz = size(z,2);

[x, y, z]   = meshgrid(x,y,z);
x           = x(:);
y           = y(:);
z           = z(:);
numnp       = length(x);

figure(1)
plot3(x, y, z, 'k.')
hold on

%%% FIND CONNECTIVITY
c   = 1;
for k = 1:1:noz-1
    for j = 1:1:nox-1
        for i = 1:1:noy-1
            n1  = i     + noy*(j-1) + nox*noy*(k-1);
            n2  = i     + noy*j     + nox*noy*(k-1);
            n3  = n2 + 1;
            n4  = n1 + 1;
            n5  = i     + noy*(j-1) + nox*noy*k;
            n6  = i     + noy*j     + nox*noy*k;
            n7  = n6 + 1;
            n8  = n5 + 1;
            
            np(c,1:8)   = [n1 n2 n3 n4 n5 n6 n7 n8];
            
            plot3(x(n1), y(n1), z(n1), 'rs')
            plot3(x(n2), y(n2), z(n2), 'gs')
            plot3(x(n3), y(n3), z(n3), 'bs')
            plot3(x(n4), y(n4), z(n4), 'ms')
            plot3(x(n5), y(n5), z(n5), 'r^')
            plot3(x(n6), y(n6), z(n6), 'g^')
            plot3(x(n7), y(n7), z(n7), 'b^')
            plot3(x(n8), y(n8), z(n8), 'm^')
            
            c   = c + 1;
        end
    end
end

% APPLY BCS
% TRACTION FREE : nps numels
% SHEAR FREE    : npsyms, numelsyms
nps     = [];
npsyms  = [];

%%% TOP FACE - TRACTION
figure(2)
subplot(2,3,6)
plot3(x, y, z, 'k.')
hold on
for k = noz-1
    for j = 1:1:nox-1
        for i = 1:1:noy-1
            n1  = i     + noy*(j-1) + nox*noy*(k-1);
            n2  = i     + noy*j     + nox*noy*(k-1);
            n3  = n2 + 1;
            n4  = n1 + 1;
            n5  = i     + noy*(j-1) + nox*noy*k;
            n6  = i     + noy*j     + nox*noy*k;
            n7  = n6 + 1;
            n8  = n5 + 1;
            
            % nps = [nps; n5 n6 n7 n8];    %%% OK
            
            plot3(x(n5), y(n5), z(n5), 'g^')
            plot3(x(n6), y(n6), z(n6), 'gv')
            plot3(x(n7), y(n7), z(n7), 'g<')
            plot3(x(n8), y(n8), z(n8), 'g>')
        end
    end
end

%%% BOTTOM FACE - TRACTION
figure(2)
subplot(2,3,1)
plot3(x, y, z, 'k.')
hold on
for k = 1
    for j = 1:1:nox-1
        for i = 1:1:noy-1
            n1  = i     + noy*(j-1) + nox*noy*(k-1);
            n2  = i     + noy*j     + nox*noy*(k-1);
            n3  = n2 + 1;
            n4  = n1 + 1;
            n5  = i     + noy*(j-1) + nox*noy*k;
            n6  = i     + noy*j     + nox*noy*k;
            n7  = n6 + 1;
            n8  = n5 + 1;
            
            % nps = [nps; n1 n4 n3 n2];    %%% OK
            
            plot3(x(n1), y(n1), z(n1), 'g^')
            plot3(x(n4), y(n4), z(n4), 'gv')
            plot3(x(n3), y(n3), z(n3), 'g<')
            plot3(x(n2), y(n2), z(n2), 'g>')
        end
    end
end

%%% FOUR SIDES - TRACTION FREE
%%% FACE 1
figure(2)
subplot(2,3,2)
plot3(x, y, z, 'k.')
hold on
for k = 1:1:noz-1
    for j = 1:1:nox-1
        for i = noy-1
            n1  = i     + noy*(j-1) + nox*noy*(k-1);
            n2  = i     + noy*j     + nox*noy*(k-1);
            n3  = n2 + 1;
            n4  = n1 + 1;
            n5  = i     + noy*(j-1) + nox*noy*k;
            n6  = i     + noy*j     + nox*noy*k;
            n7  = n6 + 1;
            n8  = n5 + 1;
            
            nps = [nps; n4 n8 n7 n3];    %%% OK
            
            plot3(x(n4), y(n4), z(n4), 'r^')
            plot3(x(n8), y(n8), z(n8), 'rv')
            plot3(x(n7), y(n7), z(n7), 'r<')
            plot3(x(n3), y(n3), z(n3), 'r>')
        end
    end
end

%%% FACE 2
subplot(2,3,3)
plot3(x, y, z, 'k.')
hold on
for k = 1:1:noz-1
    for j = 1:1:nox-1
        for i = 1
            n1  = i     + noy*(j-1) + nox*noy*(k-1);
            n2  = i     + noy*j     + nox*noy*(k-1);
            n3  = n2 + 1;
            n4  = n1 + 1;
            n5  = i     + noy*(j-1) + nox*noy*k;
            n6  = i     + noy*j     + nox*noy*k;
            n7  = n6 + 1;
            n8  = n5 + 1;
            
            nps = [nps; n1 n2 n6 n5];    %%% OK
            
            plot3(x(n1), y(n1), z(n1), 'r^')
            plot3(x(n2), y(n2), z(n2), 'rv')
            plot3(x(n6), y(n6), z(n6), 'r<')
            plot3(x(n5), y(n5), z(n5), 'r>')
        end
    end
end

%%% FACE 3
subplot(2,3,4)
plot3(x, y, z, 'k.')
hold on
for k = 1:1:noz-1
    for j = 1
        for i = 1:1:noy-1
            n1  = i     + noy*(j-1) + nox*noy*(k-1);
            n2  = i     + noy*j     + nox*noy*(k-1);
            n3  = n2 + 1;
            n4  = n1 + 1;
            n5  = i     + noy*(j-1) + nox*noy*k;
            n6  = i     + noy*j     + nox*noy*k;
            n7  = n6 + 1;
            n8  = n5 + 1;
            
            nps = [nps; n1 n5 n8 n4];    %%% OK
            
            plot3(x(n1), y(n1), z(n1), 'r^')
            plot3(x(n5), y(n5), z(n5), 'rv')
            plot3(x(n8), y(n8), z(n8), 'r<')
            plot3(x(n4), y(n4), z(n4), 'r>')
        end
    end
end

%%% FACE 4
subplot(2,3,5)
plot3(x, y, z, 'k.')
hold on
for k = 1:1:noz-1
    for j = nox-1
        for i = 1:1:noy-1
            n1  = i     + noy*(j-1) + nox*noy*(k-1);
            n2  = i     + noy*j     + nox*noy*(k-1);
            n3  = n2 + 1;
            n4  = n1 + 1;
            n5  = i     + noy*(j-1) + nox*noy*k;
            n6  = i     + noy*j     + nox*noy*k;
            n7  = n6 + 1;
            n8  = n5 + 1;
            
            nps = [nps; n2 n3 n7 n6];    %%% OK
            
            plot3(x(n2), y(n2), z(n2), 'r^')
            plot3(x(n3), y(n3), z(n3), 'rv')
            plot3(x(n7), y(n7), z(n7), 'r<')
            plot3(x(n6), y(n6), z(n6), 'r>')
        end
    end
end

numel       = size(np,1);
numels      = size(nps,1);
numelsyms   = size(npsyms,1);

nnps    = 4;
meltyp  = 8;
numeq   = 3*numnp;
nnpe    = meltyp;

% Save the meshdata
PFNAME_MESH = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_MESH);
save(PFNAME_MESH, ...
    'numnp', 'numeq', 'numel', 'numels', 'numelsyms', 'meltyp', 'nnpe', 'nnps', 'np', 'nps', 'npsyms', 'x', 'y', 'z')

figure(3)
hold on
for i = 1:1:numel
    np(i,:);
    plot3(x(np(i,1)), y(np(i,1)), z(np(i,1)), 'k.')
    plot3(x(np(i,2)), y(np(i,2)), z(np(i,2)), 'k.')
    plot3(x(np(i,3)), y(np(i,3)), z(np(i,3)), 'k.')
    plot3(x(np(i,4)), y(np(i,4)), z(np(i,4)), 'k.')
    plot3(x(np(i,5)), y(np(i,5)), z(np(i,5)), 'k.')
    plot3(x(np(i,6)), y(np(i,6)), z(np(i,6)), 'k.')
    plot3(x(np(i,7)), y(np(i,7)), z(np(i,7)), 'k.')
    plot3(x(np(i,8)), y(np(i,8)), z(np(i,8)), 'k.')
    plot3(mean(x(np(i,:))), mean(y(np(i,:))), mean(z(np(i,:))), 'y.')
end

for i = 1:1:numelsyms
    npsyms(i,:);
    for j = 1:1:4
        plot3(x(npsyms(i,j)), y(npsyms(i,j)), z(npsyms(i,j)), 'ro')
    end
end

for i = 1:1:numels
    for j = 1:1:4
        plot3(x(nps(i,j)), y(nps(i,j)), z(nps(i,j)), 'gs')
        % pause
    end
end

for i = 1:1:InputData.ndv
    plot3(...
        InputData.dvx(i), ...
        InputData.dvy(i), ...
        InputData.dvz(i), ...
        'b^')
    % pause(0.5)
end
disp('done')