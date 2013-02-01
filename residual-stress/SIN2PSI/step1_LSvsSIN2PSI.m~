clear all
close all
clc

pname_SPF   = '/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode3.Ti811/Final_SPF';

angle_pos   = [0 90]';
radial_pos  = 1:1:21;

na  = length(angle_pos);
nr  = length(radial_pos);
LSRange = [0 0];
for i = 1:1:na
    for j = 1:1:nr
        fname_SPF   = ['Alp', num2str(angle_pos(i)), ...
            '_pos', num2str(radial_pos(j)), ...
            'LP2smooth.res.mat'];
        pfname  = fullfile(pname_SPF, fname_SPF);
        SPFData = load(pfname);
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%% 110
        q   = SPFData.r110(:,1:3);
        LS  = SPFData.r110(:,4);
        
        if angle_pos(i) == 0
            idx = find(q(:,2) == 0);
        elseif angle_pos(i) == 90
            idx = find(q(:,2) == 0);
        else
            disp('omg angle not considered')
        end
        LS110{i,j}  = [q(idx,:) LS(idx,:)];
        
        if min(LS(idx,:)) < LSRange(1)
            LSRange(1) = min(LS(idx,:));
        end
        if max(LS(idx,:)) > LSRange(2)
            LSRange(2)  = max(LS(idx,:));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%% 112
        q   = SPFData.r112(:,1:3);
        LS  = SPFData.r112(:,4);
        
        if angle_pos(i) == 0
            idx = find(q(:,2) == 0);
        elseif angle_pos(i) == 90
            idx = find(q(:,2) == 0);
        else
            disp('omg angle not considered')
        end
        LS112{i,j}  = [q(idx,:) LS(idx,:)];
        
        if min(LS(idx,:)) < LSRange(1)
            LSRange(1) = min(LS(idx,:));
        end
        if max(LS(idx,:)) > LSRange(2)
            LSRange(2)  = max(LS(idx,:));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT SIN2PSI VS. STRAIN
for i = 1:1:na
    for j = 1:1:nr
        %%%%%%%%%%%%%%%%%%%%%%
        %%% 110
        q   = LS110{i,j}(:,1:3);
        LS  = LS110{i,j}(:,4);
        sin2psi = 1 - q(:,3).^2;

        figure(j)
        subplot(2,2,1+(i-1)*2)
        plot(sin2psi, LS, 'k.')
        axis([0 1 LSRange])
        title(['\{11.0\}'])
        grid on

        %%%%%%%%%%%%%%%%%%%%%%
        %%% 112
        q   = LS112{i,j}(:,1:3);
        LS  = LS112{i,j}(:,4);
        sin2psi = 1 - q(:,3).^2;

        figure(j)
        subplot(2,2,2+(i-1)*2)
        plot(sin2psi, LS, 'k.')
        axis([0 1 LSRange])
        title(['\{11.2\}'])
        grid on
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE SIGMAXX, SIGMAYY FROM SIN2PSI
% ndiv    = 1000;
% qsym    = HexSymmetries;
% 
% E_hkl   = 120e3;    % MPa (BULK)
% nu_hkl  = 0.32;
% 
% S1  = -nu_hkl/E_hkl;
% S2  = 2*(1 + nu_hkl)/E_hkl;
% 
% p   = [ ...
%     1 0 0; ...
%     ]';
% for j = 1:1:nr
%     for i = 1:1:na       
%         %%%%%%%%%%%%%%%%%%%%%%
%         %%% 100
%         q   = LS100{i,j}(:,1:3);
%         LS  = LS100{i,j}(:,4);
%         sin2psi = 1 - q(:,3).^2;
%         
%         
%         p100{i,j}       = polyfit(sin2psi, LS, 1);
%         p100ft{i,j}     = polyfit(sin2psi(1:4:end), LS(1:4:end), 1);
%         
%         idx = sin2psi < 0.6;
%         p100ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);
%         
%         sx100(i,j)      = p100{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy100(i,j)      = -p100{i,j}(2)*E_hkl/nu_hkl-sx100(i,j);
%         sxft100(i,j)    = p100ft{i,j}(1)*E_hkl/(1 + nu_hkl);
%         syft100(i,j)    = -p100ft{i,j}(2)*E_hkl/nu_hkl-sxft100(i,j);
%         sx100ft2(i,j)   = p100ft2{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy100ft2(i,j)   = -p100ft2{i,j}(2)*E_hkl/nu_hkl-sx100ft2(i,j);
%         
%         %%%%%%%%%%%%%%%%%%%%%%
%         %%% 002
%         q   = LS002{i,j}(:,1:3);
%         LS  = LS002{i,j}(:,4);
%         sin2psi = 1 - q(:,3).^2;
%         
%         p002{i,j}   = polyfit(sin2psi, LS, 1);
%         p002ft{i,j} = polyfit(sin2psi(1:4:end), LS(1:4:end), 1);
%         
%         idx = sin2psi < 0.6;
%         p002ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);
%         
%         sx002(i,j)      = p002{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy002(i,j)      = -p002{i,j}(2)*E_hkl/nu_hkl-sx002(i,j);
%         sxft002(i,j)    = p002ft{i,j}(1)*E_hkl/(1 + nu_hkl);
%         syft002(i,j)    = -p002ft{i,j}(2)*E_hkl/nu_hkl-sxft002(i,j);
%         sx002ft2(i,j)   = p002ft2{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy002ft2(i,j)   = -p002ft2{i,j}(2)*E_hkl/nu_hkl-sx002ft2(i,j);
%         
%         %%%%%%%%%%%%%%%%%%%%%%
%         %%% 101
%         q   = LS101{i,j}(:,1:3);
%         LS  = LS101{i,j}(:,4);
%         sin2psi = 1 - q(:,3).^2;
%         
%         p101{i,j}   = polyfit(sin2psi, LS, 1);
%         p101ft{i,j} = polyfit(sin2psi(1:4:end), LS(1:4:end), 1);
%         
%         idx = sin2psi < 0.6;
%         p101ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);
%         
%         sx101(i,j)      = p101{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy101(i,j)      = -p101{i,j}(2)*E_hkl/nu_hkl-sx101(i,j);
%         sxft101(i,j)    = p101ft{i,j}(1)*E_hkl/(1 + nu_hkl);
%         syft101(i,j)    = -p101ft{i,j}(2)*E_hkl/nu_hkl-sxft101(i,j);
%         sx101ft2(i,j)   = p101ft2{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy101ft2(i,j)   = -p101ft2{i,j}(2)*E_hkl/nu_hkl-sx101ft2(i,j);
%         
%         %%%%%%%%%%%%%%%%%%%%%%
%         %%% 110
%         q   = LS110{i,j}(:,1:3);
%         LS  = LS110{i,j}(:,4);
%         sin2psi = 1 - q(:,3).^2;
%         
%         p110{i,j}   = polyfit(sin2psi, LS, 1);
%         p110ft{i,j} = polyfit(sin2psi(1:4:end), LS(1:4:end), 1);
%         
%         idx = sin2psi < 0.6;
%         p110ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);
%         
%         sx110(i,j)      = p110{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy110(i,j)      = -p110{i,j}(2)*E_hkl/nu_hkl-sx110(i,j);
%         sxft110(i,j)    = p110ft{i,j}(1)*E_hkl/(1 + nu_hkl);
%         syft110(i,j)    = -p110ft{i,j}(2)*E_hkl/nu_hkl-sxft110(i,j);
%         sx110ft2(i,j)   = p110ft2{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy110ft2(i,j)   = -p110ft2{i,j}(2)*E_hkl/nu_hkl-sx110ft2(i,j);
%         
%         %%%%%%%%%%%%%%%%%%%%%%
%         %%% 112
%         q   = LS112{i,j}(:,1:3);
%         LS  = LS112{i,j}(:,4);
%         sin2psi = 1 - q(:,3).^2;
%         
%         p112{i,j}   = polyfit(sin2psi, LS, 1);
%         p112ft{i,j} = polyfit(sin2psi(1:4:end), LS(1:4:end), 1);
%         
%         idx = sin2psi < 0.6;
%         p112ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);
%         
%         sx112(i,j)      = p112{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy112(i,j)      = -p112{i,j}(2)*E_hkl/nu_hkl-sx112(i,j);
%         sxft112(i,j)    = p112ft{i,j}(1)*E_hkl/(1 + nu_hkl);
%         syft112(i,j)    = -p112ft{i,j}(2)*E_hkl/nu_hkl-sxft112(i,j);
%         sx112ft2(i,j)   = p112ft2{i,j}(1)*E_hkl/(1 + nu_hkl);
%         sy112ft2(i,j)   = -p112ft2{i,j}(2)*E_hkl/nu_hkl-sx112ft2(i,j);
%     end
% end
% 
% DataRange(1)    = min([ ...
%     sx100(:); sxft100(:); sx100ft2(:); ...
%     sx002(:); sxft002(:); sx002ft2(:); ...
%     sx101(:); sxft101(:); sx101ft2(:); ...
%     sx110(:); sxft110(:); sx110ft2(:); ...
%     sx112(:); sxft112(:); sx112ft2(:); ...
%     sy100(:); syft100(:); sy100ft2(:); ...
%     sy002(:); syft002(:); sy002ft2(:); ...
%     sy101(:); syft101(:); sy101ft2(:); ...
%     sy110(:); syft110(:); sy110ft2(:); ...
%     sy112(:); syft112(:); sy112ft2(:); ...
%     ]);
% DataRange(2)    = max([ ...
%     sx100(:); sxft100(:); sx100ft2(:); ...
%     sx002(:); sxft002(:); sx002ft2(:); ...
%     sx101(:); sxft101(:); sx101ft2(:); ...
%     sx110(:); sxft110(:); sx110ft2(:); ...
%     sx112(:); sxft112(:); sx112ft2(:); ...
%     sy100(:); syft100(:); sy100ft2(:); ...
%     sy002(:); syft002(:); sy002ft2(:); ...
%     sy101(:); syft101(:); sy101ft2(:); ...
%     sy110(:); syft110(:); sy110ft2(:); ...
%     sy112(:); syft112(:); sy112ft2(:); ...
%     ]);
% 
% SOL     = load('/home/jspark/Documents/work/prjTi-RS/Atish/Stress_code-2D/STRESS.SOL.mat');
% X0      = SOL.X(:,17);
% Y0      = SOL.Y(:,17);
% R       = sqrt(X0.^2 + Y0.^2);
% SIGXX   = [ ...
%     SOL.Dx(:,17)'; ...
%     SOL.Dx(:,9)'; ...
%     ];
% SIGYY   = [ ...
%     SOL.Dy(:,17)'; ...
%     SOL.Dy(:,9)'; ...
%     ];
% 
% x   = [6.5000 6.6350 6.7770 6.9300 7.0950 7.2700 7.4600 7.6650 7.8900 8.1300 ...
%     8.4100 8.7100 9.0500 9.4200 9.8500 10.3500 10.9500 11.6500 12.5500 ...
%     13.6000 15.0500];
% for i = 1:1:2
%     fig = figure(1);
%     orient rotated
%     set(fig, 'position', [1681 27 1280 940], ...
%         'PaperOrientation', 'portrait')
%     
%     subplot(2,5,1+5*(i-1))
%     plot(x,sx100(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'b')
%     hold on
%     plot(x,sxft100(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'r')
%     plot(x,sx100ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'g')
%     plot(x,sy100(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'none')
%     plot(x,syft100(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'none')
%     plot(x,sy100ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'none')
%     plot(R,SIGXX(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'k')
%     plot(R,SIGYY(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'none')
%     axis([min(x) max(x) DataRange])
%     title(['\{10.0\}'])
%     xlabel('radial position (mm)')
%     ylabel('stress (MPa)')
%     grid on
%     
%     subplot(2,5,2+5*(i-1))
%     plot(x,sx002(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'b')
%     hold on
%     plot(x,sxft002(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'r')
%     plot(x,sx002ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'g')
%     plot(x,sy002(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'none')
%     plot(x,syft002(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'none')
%     plot(x,sy002ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'none')
%     plot(R,SIGXX(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'k')
%     plot(R,SIGYY(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'none')
%     axis([min(x) max(x) DataRange])
%     title(['\{00.2\}'])
%     xlabel('radial position (mm)')
%     ylabel('stress (MPa)')
%     grid on
%     
%     subplot(2,5,3+5*(i-1))
%     plot(x,sx101(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'b')
%     hold on
%     plot(x,sxft101(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'r')
%     plot(x,sx101ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'g')
%     plot(x,sy101(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'none')
%     plot(x,syft101(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'none')
%     plot(x,sy101ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'none')
%     plot(R,SIGXX(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'k')
%     plot(R,SIGYY(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'none')
%     axis([min(x) max(x) DataRange])
%     title(['\{10.1\}'])
%     xlabel('radial position (mm)')
%     ylabel('stress (MPa)')
%     grid on
%     
%     subplot(2,5,4+5*(i-1))
%     plot(x,sx110(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'b')
%     hold on
%     plot(x,sxft110(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'r')
%     plot(x,sx110ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'g')
%     plot(x,sy110(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'none')
%     plot(x,syft110(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'none')
%     plot(x,sy110ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'none')
%     plot(R,SIGXX(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'k')
%     plot(R,SIGYY(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'none')
%     axis([min(x) max(x) DataRange])
%     title(['\{11.0\}'])
%     xlabel('radial position (mm)')
%     ylabel('stress (MPa)')
%     grid on
%     
%     subplot(2,5,5+5*(i-1))
%     plot(x,sx112(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'b')
%     hold on
%     plot(x,sxft112(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'r')
%     plot(x,sx112ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'g')
%     plot(x,sy112(i,:), 'o', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'b', ...
%         'MarkerFaceColor', 'none')
%     plot(x,syft112(i,:), '^', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'none')
%     plot(x,sy112ft2(i,:), 'v', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'g', ...
%         'MarkerFaceColor', 'none')
%     plot(R,SIGXX(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'k')
%     plot(R,SIGYY(i,:), 's', ...
%         'MarkerSize', 8, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', 'none')
%     axis([min(x) max(x) DataRange])
%     title(['\{11.2\}'])
%     xlabel('radial position (mm)')
%     ylabel('stress (MPa)')
%     grid on
% end
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPTIMIZE E / nu
SOL     = load('/home/jspark/Documents/work/prjTi-RS/Atish/Stress_code-2D/STRESS.SOL.mat');
X0      = SOL.X(:,17);
Y0      = SOL.Y(:,17);
R       = sqrt(X0.^2 + Y0.^2);
SIGXX   = [ ...
    SOL.Dx(:,17)'; ...
    SOL.Dx(:,9)'; ...
    ];
SIGXX   = (SIGXX(:,1:end-1) + SIGXX(:,2:end))./2;

SIGYY   = [ ...
    SOL.Dy(:,17)'; ...
    SOL.Dy(:,9)'; ...
    ];
SIGYY   = (SIGYY(:,1:end-1) + SIGYY(:,2:end))./2;

SIGXY   = [ ...
    SOL.Dxy(:,17)'; ...
    SOL.Dxy(:,9)'; ...
    ];
SIGXY   = (SIGXY(:,1:end-1) + SIGXY(:,2:end))./2;

x   = [6.5000 6.6350 6.7770 6.9300 7.0950 7.2700 7.4600 7.6650 7.8900 8.1300 ...
    8.4100 8.7100 9.0500 9.4200 9.8500 10.3500 10.9500 11.6500 12.5500 ...
    13.6000 15.0500];

for j = 1:1:nr
    for i = 1:1:na       
        %%%%%%%%%%%%%%%%%%%%%%
        %%% hkl
        q   = LS100{i,j}(:,1:3);
        LS  = LS100{i,j}(:,4);
        
%         q   = LS002{i,j}(:,1:3);
%         LS  = LS002{i,j}(:,4);
        
%         q   = LS101{i,j}(:,1:3);
%         LS  = LS101{i,j}(:,4);
        
%         q   = LS110{i,j}(:,1:3);
%         LS  = LS110{i,j}(:,4);
        
%         q   = LS112{i,j}(:,1:3);
%         LS  = LS112{i,j}(:,4);
        
        sin2psi = 1 - q(:,3).^2;
        
%         p100{i,j}       = polyfit(sin2psi, LS, 1);
%         m100(i,j)       = p100{i,j}(1);
%         b100(i,j)       = p100{i,j}(2);
        
        idx = sin2psi < 0.6;
        p100ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);       
        m100(i,j)       = p100ft2{i,j}(1);
        b100(i,j)       = p100ft2{i,j}(2);
    end
end

%%% ALL DATA
xdata   = [m100(:) b100(:)];
ydata   = [SIGXX(:); SIGYY(:)]./1000;

% %%% SIGXX - A0
xdata   = [m100(1,:); b100(1,:)]';
ydata   = [SIGXX(1,:)]'./1000;

%%% SIGXX - A90
xdata   = [m100(2,:); b100(2,:)]';
ydata   = [SIGXX(2,:)]'./1000;

%%% SIGYY - A0
xdata   = [m100(1,:); b100(1,:)]';
ydata   = [SIGYY(1,:)]'./1000;

%%% SIGYY - A90
xdata   = [m100(2,:); b100(2,:)]';
ydata   = [SIGYY(2,:)]'./1000;

E12_100i    = [100 0.3]';

%%% NEED TO CHANGE "OptimizeProperty" DEPENDING ON CASES
yfunc0      = OptimizeProperty(E12_100i,xdata);

LB      = [-inf 0.17];
UB      = [inf 0.47];
[E12_100f,~,~,ef,output,lambda,~]   = lsqcurvefit(@OptimizeProperty, E12_100i, xdata, ydata, ...
    LB, UB);
yfunc   = OptimizeProperty(E12_100f,xdata);

E12_100f(1)
E12_100f(2)

figure,
plot(ydata, 'k.')
hold on
plot(yfunc0, 'b.')
plot(yfunc, 'r.')
