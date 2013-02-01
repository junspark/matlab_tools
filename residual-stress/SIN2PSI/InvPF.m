clc
clear all
close all

% Input r_E value and Poisson's ratio
r_E = 0;
nu  = 0;

% Flags to plot figures
plotEd = 0;
plotEd_hkl = 0;
plotTaylor = 0;
plotSchmid = 0;
plotInvSchmid = 0;
plotS2S = 0;

%__________________________________________________________________________
%

if nu == 0.35
    if r_E == 1.0
        C11 = 110.700; C12 = 59.200; C44 = 51.500/2;
    elseif r_E == 1.2
        C11 = 107.3;   C12 = 60.9;   C44 = 56.6/2;
    elseif r_E == 1.3
        C11 = 105.729; C12 = 61.686; C44 = 58.957/2;
    elseif r_E == 1.4
        C11 = 104.356; C12 = 62.372; C44 = 61.015/2;
    elseif r_E == 1.5
        C11 = 103.099; C12 = 63.00;  C44 = 62.901/2;
    elseif r_E == 1.7
        C11 = 100.88;  C12 = 64.11;  C44 = 66.23/2;
    elseif r_E == 2.0
        C11 = 98.141;  C12 = 65.479; C44 = 70.338/2;
    elseif r_E == 2.5
        C11 = 94.691;  C12 = 67.204; C44 = 75.513/2;
    elseif r_E == 3.0
        C11 = 92.166;  C12 = 68.467; C44 = 79.301/2;
    end
elseif nu == 0.25
    if r_E == 1.2
        C11 = 79.381;  C12 = 29.754; C44 = 61.485/2;
    end
else
    C11 = 166;  C12 = 120; C44 = 76;
end

%__________________________________________________________________________
%

% Number of contour lines
ncont = 20;

% Stiffness matrix (6x6)
C = [C11 C12 C12 0   0   0
    C12 C11 C12 0   0   0
    C12 C12 C11 0   0   0
    0   0   0   C44 0   0
    0   0   0   0   C44 0
    0   0   0   0   0   C44];

% Compliance matrix (6x6)
S=inv(C);

% Applied stress in sample frame
Stress_xyz = [100 0 0;
    0 0 0;
    0 0 0];

% Generate grid of theta and phi points
% tan(phi) = (z/y)*sin(theta)
% z/y = 1 (why?)
% phi is angle with x-y plane, not angle with vertical
thetaGrid   = [0:0.3:45];
thetaGrid   = thetaGrid*pi/180;
phiGrid     = atan(sin(thetaGrid));
[theta, phi]    = meshgrid(thetaGrid, phiGrid);
theta   = reshape(theta, 1, length(thetaGrid)^2);
phi     = reshape(phi, 1, length(thetaGrid)^2);
eta     = 0;

m = cos(phi).*(cos(theta) + sin(theta)) - sin(phi);
m = (cos(phi).*cos(theta) + sin(phi)).*m;
m = m/sqrt(6);

Ed      = 0*theta;
taylor  = 0*theta;

% rotation about x'' by eta  x''y''z'' => 123
Rh  = RMatOfQuat(QuatOfAngleAxis(eta*pi/180, [1 0 0]'));          % identity matrix

R = zeros(3,3,length(theta));

% xstal stays in xyz (xyz is sample system)
% loading is moved to 123 (123 is crystal system)
for i = 1:1:length(theta)

    % Rotation matrix that takes sample coordinate system to crystal
    % coordinate system
    
    % rotate about y by phi xyz => x'y'z'
    Rf  = [
        cos(phi(i)) 0 sin(phi(i));
        0 1 0;
        -sin(phi(i)) 0 cos(phi(i));
        ];

    % rotation about z' by -theta x'y'z' => x''y''z''
    Rq  = [
        cos(-theta(i)) -sin(-theta(i)) 0;
        sin(-theta(i)) cos(-theta(i)) 0;
        0 0 1;
        ];

    % R: xyz => 123 (sample to crystal)
    R(:,:,i) = Rh*Rq*Rf;

    % transform from 123 => xyz
    % stress is in crystal system
    %Stress_123  = R(:,:,i)'*Stress_xyz*R(:,:,i);
    Stress_123  = R(:,:,i)*Stress_xyz*R(:,:,i)';

    Stress_123  = [
        Stress_123(1,1);
        Stress_123(2,2);
        Stress_123(3,3);
        Stress_123(2,3);
        Stress_123(1,3);
        Stress_123(1,2);
        ];

    % Strain in crystal system
    Strain_123   = S*Stress_123;
    Strain_123   = [
        Strain_123(1,1) Strain_123(6,1)/2 Strain_123(5,1)/2;
        Strain_123(6,1)/2 Strain_123(2,1) Strain_123(4,1)/2;
        Strain_123(5,1)/2 Strain_123(4,1)/2 Strain_123(3,1);
        ];

    % transform from xyz => 123
    % rotate strain in crystal system to sample system
    %Strain_xyz      = R(:,:,i)*Strain_123*R(:,:,i)';
    Strain_xyz      = R(:,:,i)'*Strain_123*R(:,:,i);

    Ed(i)   = Stress_xyz(1,1)/Strain_xyz(1,1);  % stiffness in GPa

    poissonxy(i) = -Strain_xyz(2,2)/Strain_xyz(1,1);
    poissonxz(i) = -Strain_xyz(3,3)/Strain_xyz(1,1);
end
return
%__________________________________________________________________________
%

load taylor_m02_R03.mat

i       = (tan(phi) > sin(theta));
Ed(i)   = NaN;
taylor(i) = NaN;
m(i) = NaN;
poissonxy(i) = NaN;
poissonxz(i) = NaN;

theta   = reshape(theta, length(thetaGrid), length(thetaGrid));
phi     = reshape(phi, length(thetaGrid), length(thetaGrid));
Ed      = reshape(Ed, length(thetaGrid), length(thetaGrid));
taylor = reshape(taylor, length(thetaGrid), length(thetaGrid));
m      = reshape(m, length(thetaGrid), length(thetaGrid));
poissonxy = reshape(poissonxy, length(thetaGrid), length(thetaGrid));
poissonxz = reshape(poissonxz, length(thetaGrid), length(thetaGrid));


figure('Color',[1 1 1])
hold on
contourf(theta*180/pi, phi*180/pi, poissonxy, ncont)
plot(thetaGrid*180/pi, atan(sin(thetaGrid))*180/pi,'Color','k','LineStyle','-','LineWidth',2)
axis tight; axis square
set(gca,'FontSize',22,'FontName','Times','FontWeight','Normal')
xlabel('   (\circ)')
ylabel('          (\circ)','Rotation',360)
%text('Interpreter','latex','String','$\theta$','Position',[20.2 -3.7],'FontSize',24)
%text('Interpreter','latex','String','$\phi$','Position',[48.5 16.6],'FontSize',24)
set(colorbar,'FontSize',22,'FontName','Times','FontWeight','Normal','Location','WestOutside')
set(gca,'Box','off','LineWidth',2,'FontName','Times','FontSize',22,'FontWeight','Normal','YAxisLocation','Right')
%text(-11.5,35,'GPa','FontSize',22,'FontName','Times','FontWeight','Normal')


figure('Color',[1 1 1])
hold on
contourf(theta*180/pi, phi*180/pi, poissonxz, ncont)
plot(thetaGrid*180/pi, atan(sin(thetaGrid))*180/pi,'Color','k','LineStyle','-','LineWidth',2)
axis tight; axis square
set(gca,'FontSize',22,'FontName','Times','FontWeight','Normal')
xlabel('   (\circ)')
ylabel('          (\circ)','Rotation',360)
%text('Interpreter','latex','String','$\theta$','Position',[20.2 -3.7],'FontSize',24)
%text('Interpreter','latex','String','$\phi$','Position',[48.5 16.6],'FontSize',24)
set(colorbar,'FontSize',22,'FontName','Times','FontWeight','Normal','Location','WestOutside')
set(gca,'Box','off','LineWidth',2,'FontName','Times','FontSize',22,'FontWeight','Normal','YAxisLocation','Right')
%text(-11.5,35,'GPa','FontSize',22,'FontName','Times','FontWeight','Normal')


return


%--------------------------------------------------------------------------
%                         Directonal modulus (Ed)
%--------------------------------------------------------------------------

if plotEd
    figure('Color',[1 1 1])
    hold on
    contourf(theta*180/pi, phi*180/pi, Ed, ncont)
    plot(thetaGrid*180/pi, atan(sin(thetaGrid))*180/pi,'Color','k','LineStyle','-','LineWidth',2)
    axis tight; axis square
    set(gca,'FontSize',22,'FontName','Times','FontWeight','Normal')
    xlabel('   (\circ)')
    ylabel('          (\circ)','Rotation',360)
    text('Interpreter','latex','String','$\theta$','Position',[20.2 -3.7],'FontSize',24)
    text('Interpreter','latex','String','$\phi$','Position',[48.5 16.6],'FontSize',24)
    set(colorbar,'FontSize',22,'FontName','Times','FontWeight','Normal','Location','WestOutside')
    set(gca,'Box','off','LineWidth',2,'FontName','Times','FontSize',22,'FontWeight','Normal','YAxisLocation','Right')
    text(-11.5,35,'GPa','FontSize',22,'FontName','Times','FontWeight','Normal') 
    % saveas(gcf,'~/Desktop/test.eps','psc2')
end

%--------------------------------------------------------------------------

% Plot directonal modulus with hkl
if plotEd_hkl
    figure; hold on
    set(gcf,'Color',[1 1 1])
    caxis([33 102])
    colorbar('FontSize', 16, 'FontWeight', 'Bold', 'Location','WestOutside','CLim',[33 102])
    contourf(theta*180/pi, phi*180/pi, Ed, ncont)
    plot(thetaGrid*180/pi, atan(sin(thetaGrid))*180/pi, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2)
    set(gca, 'Box', 'off', 'LineWidth', 2, 'FontSize', 16, 'FontWeight', 'bold', 'YAxisLocation', 'right')
    text(-17,36,'GPa','FontSize',16,'FontWeight','Bold')
    plot(0:1:45, zeros(length(0:1:45)), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4)
    plot(45*ones(length(0:1:35)), 0:1:35, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4)
    title(['r_E = ', num2str(r_E,'%2.1f')], 'FontSize', 16, 'FontWeight', 'Bold')
    axis([0 45 0 37])
    axis square
    axis off

    h = 1; k = 0; l = 0;
    [theta1,phi1,r] = cart2sph(h,k,l);
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',35)
    text(-3,-2,'[100]','FontSize',16,'FontWeight','Bold')
    Ed = Ed_Eval([h k l], r_E)
    
    h = 1; k = 1; l = 0;
    [theta1,phi1,r] = cart2sph(h,k,l);
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',35)
    text(42,-2,'[110]','FontSize',16,'FontWeight','Bold')
    Ed = Ed_Eval([h k l], r_E)
    
    h = 1; k = 1; l = 1;
    [theta1,phi1,r] = cart2sph(h,k,l);
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',35)
    text(42,37,'[111]','FontSize',16,'FontWeight','Bold')
    Ed = Ed_Eval([h k l], r_E)
    
    h = 3; k = 1; l = 1;
    [theta1,phi1,r] = cart2sph(h,k,l);
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',35)
    text(12,18,'[311]','FontSize',16,'FontWeight','Bold')
    Ed = Ed_Eval([h k l], r_E)

    
end

%--------------------------------------------------------------------------

% Plot Taylor factor
if plotTaylor
    figure
    contourf(theta*180/pi, phi*180/pi, taylor, ncont)
    hold on
    plot(thetaGrid*180/pi, atan(sin(thetaGrid))*180/pi, ...
        'Color', 'k', 'LineStyle', '-', 'LineWidth', 2)
    title('Taylor factor, M', 'FontSize', 16, 'FontWeight', 'bold')
    xlabel('\theta (deg)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('\phi (deg)', 'FontSize', 16, 'FontWeight', 'bold');
    set(gca, 'Box', 'off', 'LineWidth', 2, 'FontSize', 16, 'FontWeight', 'bold', 'YAxisLocation', 'right')
    set(colorbar, 'FontSize', 16, 'FontWeight', 'bold', 'Location', 'WestOutside')
    axis square
end

%--------------------------------------------------------------------------
%                               Schmid factor
%--------------------------------------------------------------------------
if plotSchmid

    figure('Color',[1 1 1])
    hold on
    contourf(theta*180/pi,phi*180/pi, m, ncont)
    plot(thetaGrid*180/pi,atan(sin(thetaGrid))*180/pi,'Color','k','LineStyle','-','LineWidth',2)
    axis tight; axis square
    set(gca,'FontSize',22,'FontName','Times','FontWeight','Normal')
    xlabel('   (\circ)')
    ylabel('          (\circ)','Rotation',360)
    text('Interpreter','latex','String','$\theta$','Position',[20.2 -3.7],'FontSize',24)
    text('Interpreter','latex','String','$\phi$','Position',[48.5 16.6],'FontSize',24)
    set(colorbar, 'FontSize', 22, 'FontWeight', 'Normal', 'Location', ...
        'WestOutside','YTick',[0.28:0.04:0.50],'YTickLabel',{'0.28','0.32','0.36','0.40','0.44','0.48'})
    set(gca,'Box','off','LineWidth',2,'FontName','Times','FontSize', 22,...
        'FontWeight','Normal','YAxisLocation', 'right')
    
end

%--------------------------------------------------------------------------

% Plot 1/Schmid factor
if plotInvSchmid

    figure('Color',[1 1 1])
    hold on
    contourf(theta*180/pi,phi*180/pi, 1./m, ncont)
    plot(thetaGrid*180/pi,atan(sin(thetaGrid))*180/pi,'Color','k','LineStyle','-','LineWidth',2)
    axis tight; axis square
    set(gca,'FontSize',22,'FontName','Times','FontWeight','Normal')
    xlabel('   (\circ)')
    ylabel('          (\circ)','Rotation',360)
    text('Interpreter','latex','String','$\theta$','Position',[20.2 -3.7],'FontSize',24)
    text('Interpreter','latex','String','$\phi$','Position',[48.5 16.6],'FontSize',24)
    set(colorbar, 'FontSize', 22, 'FontWeight', 'Normal', 'Location', ...
        'WestOutside','YTick',[0:0.2:5],'YTickLabel',{'2.0','2.2','2.4','2.6','2.8','3.0','3.2','3.4','3.6'})
    set(gca,'Box','off','LineWidth',2,'FontName','Times','FontSize', 22,...
        'FontWeight','Normal','YAxisLocation', 'right','XTick',[0:10:50])

end

%--------------------------------------------------------------------------
%                       Strength-to-stiffness (S2S)
%--------------------------------------------------------------------------
if plotS2S
    
    figure('Color',[1 1 1])
    hold on
    %colormap Gray
    caxis([0.0286 0.0725])
    colorbar('FontSize', 22, 'FontWeight', 'Bold', 'Location','WestOutside','CLim',[0.0286 0.0725],'YTick',[0.030:0.005:0.070],'YTickLabel',{'0.030','0.035','0.040','0.045','0.050','0.055','0.060','0.065','0.070'})
    contourf(theta*180/pi, phi*180/pi, 1./(m.*Ed), ncont)
    plot(thetaGrid*180/pi, atan(sin(thetaGrid))*180/pi, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2)
    plot(0:1:45, zeros(length(0:1:45)), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3)
    plot(45*ones(length(0:1:35)), 0:1:35, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3)
    axis tight; axis square; axis off;
    set(gca,'FontSize',22,'FontName','Times','FontWeight','Normal')
    xlabel('   (\circ)')
    ylabel('          (\circ)','Rotation',360)
    text('Interpreter','latex','String','$\theta (^\circ)$','Position',[20.2 -2.5],'FontSize',24)
    text('Interpreter','latex','String','$\phi$','Position',[46.5 16.6],'FontSize',24)
    text('Interpreter','latex','String','$(^\circ)$','Position',[48.5 16.6],'FontSize',24)
    text(-11,36.5,'GPa^{-1}','FontSize',22,'FontName','Times')

    hkl = [1 0 0];
    [theta1,phi1,r] = cart2sph(hkl(1),hkl(2),hkl(3));
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',40)
    text(-4,-1.5,'\langle100\rangle','FontSize',22,'FontName','Times')
    %[schmid_to_stiffness] = S2SEval(hkl, r_E);
    %text(-4,-2,sprintf('%0.4f',schmid_to_stiffness),'FontSize',16,'FontWeight','Bold')

    hkl = [1 1 1];
    [theta1,phi1,r] = cart2sph(hkl(1),hkl(2),hkl(3));
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',40)
    text(42.5,37,'\langle111\rangle','FontSize',22,'FontName','Times')
    %[schmid_to_stiffness] = S2SEval(hkl, r_E);
    %text(42,37,sprintf('%0.4f',schmid_to_stiffness),'FontSize',16,'FontWeight','Bold')
    
    hkl = [1 1 0];
    [theta1,phi1,r] = cart2sph(hkl(1),hkl(2),hkl(3));
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',40)
    text(42.5,-1.5,'\langle110\rangle','FontSize',22,'FontName','Times')
    %[schmid_to_stiffness] = S2SEval(hkl, r_E);
    %text(42,-2,sprintf('%0.4f',schmid_to_stiffness),'FontSize',16,'FontWeight','Bold')
    
    hkl = [3 1 1];
    [theta1,phi1,r] = cart2sph(hkl(1),hkl(2),hkl(3));
    plot(theta1*180/pi,phi1*180/pi,'.','Color',[0.5 0.5 0.5],'MarkerSize',40)
    text(11.5,18,'\langle311\rangle','FontSize',22,'FontName','Times')
    %[schmid_to_stiffness] = S2SEval(hkl, r_E);
    %text(10,17.5,sprintf('%0.4f',schmid_to_stiffness),'FontSize',16,'FontWeight','Bold')
    
end


disp(['E_111/E_100 = ',num2str(max(max(Ed))/min(min(Ed)))])
disp(['Max Taylor factor over Ed = ',sprintf('%0.4f',max(max(taylor./Ed)))])
disp(['Min Taylor factor over Ed = ',sprintf('%0.4f',min(min(taylor./Ed)))])
disp(['Max Schmid factor over Ed = ',sprintf('%0.5f',max(max(1./(m.*Ed))))])
disp(['Min Schmid factor over Ed = ',sprintf('%0.5f',min(min(1./(m.*Ed))))])

