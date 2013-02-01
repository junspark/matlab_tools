clear all
close all
clc

load Residual_Stress_3D_MLSEL_v1/meshdata.mat
load('/home/jspark/Documents/work/Data/APS201103/GE/A1234.N.mat')
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode3/A1234.K.Q.MLS.0,0,0,6.mat')
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode3/A1234.SH.NEWRESULTS.mat')
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode/A1234.Keq.mat')   %%%
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode/A1234.Kfs.mat')   %%%
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode/A1234.Ksym.mat')  %%%

%%%
% MACRO STRESS
nn  = length(x);
S   = SH(1:nn*6,1);

SEQ     = Keq*S; 
SEQ     = reshape(SEQ,6,nn)';

SFS     = Kfs*S;
SFS     = reshape(SFS,6,nn)';

SSYM    = Ksym*S;
SSYM    = reshape(SSYM,6,nn)';

S11 = S(1:6:nn*6-5);
S22 = S(2:6:nn*6-4);
S33 = S(3:6:nn*6-3);
S23 = S(4:6:nn*6-2);
S13 = S(5:6:nn*6-1);
S12 = S(6:6:nn*6-0);

SVM = sqrt((...
    (S11 - S22).^2 + ...
    (S22 - S33).^2 + ...
    (S33 - S11).^2 + ...
    6*S12.*S12 + ...
    6*S13.*S13 + ...
    6*S23.*S23)./2);

SEQ = sqrt((...
    (SEQ(:,1) - SEQ(:,2)).^2 + ...
    (SEQ(:,2) - SEQ(:,3)).^2 + ...
    (SEQ(:,3) - SEQ(:,1)).^2 + ...
    6*SEQ(:,4).*SEQ(:,4) + ...
    6*SEQ(:,5).*SEQ(:,5) + ...
    6*SEQ(:,6).*SEQ(:,6))./2);

disp([norm(SEQ) norm(SFS) norm(SSYM)])
disp([norm(Keq,1) norm(Kfs,1) norm(Ksym,1)])

f   = figure(1000);
c   = zeros(nn,1);
set(f, 'Position', [110 30 1570 915])
view([45 30])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on');
axis equal tight off off
caxis([-1 1])

c   = S11;%%
f   = figure(1);
set(f, 'Position', [110 30 1570 915])
subplot(2,3,1)
title('S_{11}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])

c   = S22;%%
f   = figure(1);
set(f, 'Position', [110 30 1570 915])
subplot(2,3,2)
title('S_{22}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])

c   = S33;%%
f   = figure(1);
set(f, 'Position', [110 30 1570 915])
subplot(2,3,3)
title('S_{33}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])

c   = S12;%%
f   = figure(1);
set(f, 'Position', [110 30 1570 915])
subplot(2,3,4)
title('S_{12}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])

c   = S13;%%
f   = figure(1);
set(f, 'Position', [110 30 1570 915])
subplot(2,3,5)
title('S_{13}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])

c   = S23;%%
f   = figure(1);
set(f, 'Position', [110 30 1570 915])
subplot(2,3,6)
title('S_{23}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])
% return
c   = SVM;%%
f   = figure(2);
title('S_{VM}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])

return
c   = SEQ;%%
f   = figure(3);
set(f, 'Position', [110 30 1570 915])
subplot(1,3,1)
title('S_{EQ}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
% caxis([-700 1100])

c   = SFS;%%
f   = figure(3);
set(f, 'Position', [110 30 1570 915])
subplot(1,3,2)
title('S_{FS}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
% caxis([-700 1100])

c   = SSYM;%%
f   = figure(3);
set(f, 'Position', [110 30 1570 915])
subplot(1,3,3)
title('S_{SYM}')%%
view([45 30])
% view([0 90])
hold on
PlotMesh(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
% caxis([-700 1100])


return
%%%
% MICRO STRESS
SIGMA   = SH(nn*6+1:end,1);
SIGMA   = N*SIGMA;
SIGMA11 = SIGMA(1:6:end-5);
SIGMA22 = SIGMA(2:6:end-4);
SIGMA33 = SIGMA(3:6:end-3);
SIGMA12 = SIGMA(4:6:end-2);
SIGMA13 = SIGMA(5:6:end-1);
SIGMA23 = SIGMA(6:6:end-0);

SIGMAVM = sqrt((...
    (SIGMA11 - SIGMA22).^2 + ...
    (SIGMA22 - SIGMA33).^2 + ...
    (SIGMA33 - SIGMA11).^2 + ...
    6*SIGMA12.*SIGMA12 + ...
    6*SIGMA13.*SIGMA13 + ...
    6*SIGMA23.*SIGMA23)./2);

c   = SIGMA22;%%
figure(1)
subplot(1,2,2)
title('\Sigma_{22}')%%
view([45 30])
% view([0 90])
hold on
PlotDV(f, np, x, y, z, c, 'ShowMesh', 'on')
axis equal tight off
colorbar('location', 'SouthOutside')
caxis([-700 1100])

return
%%%
% MACRO STRESS - MICRO STRESS
f   = figure(2);
set(f, 'Position', [1032 164 580 580])
title('S_{12} - \Sigma_{12}')%%
cMacro  = S12;%%
cMicro  = SIGMA12;%%
view([45 26])
view([0 90])
hold on
for i = 1:1:numel
    n1  = np(i,1);
    n2  = np(i,2);
    n3  = np(i,3);
    n4  = np(i,4);
    n5  = np(i,5);
    n6  = np(i,6);
    n7  = np(i,7);
    n8  = np(i,8);
    
    xx  = [x(n1) x(n2) x(n3) x(n4)...
        x(n5) x(n6) x(n7) x(n8)];
    yy  = [y(n1) y(n2) y(n3) y(n4)...
        y(n5) y(n6) y(n7) y(n8)];
    zz  = [z(n1) z(n2) z(n3) z(n4)...
        z(n5) z(n6) z(n7) z(n8)];  
    
    ccMacro = [cMacro(n1) cMacro(n2) cMacro(n3) cMacro(n4)...
        cMacro(n5) cMacro(n6) cMacro(n7) cMacro(n8)];
    ccMicro = cMicro(i)*ones(1,8);
    
    cc  = ccMacro - ccMicro;
    
    line([xx(1:4) xx(1)],[yy(1:4) yy(1)],[zz(1:4) zz(1)],'Color','k');
    line([xx(5:8) xx(5)],[yy(5:8) yy(5)],[zz(5:8) zz(5)],'Color','k');
    line([xx(1) xx(5)],[yy(1) yy(5)],[zz(1) zz(5)],'Color','k');
    line([xx(2) xx(6)],[yy(2) yy(6)],[zz(2) zz(6)],'Color','k');
    line([xx(3) xx(7)],[yy(3) yy(7)],[zz(3) zz(7)],'Color','k');
    line([xx(4) xx(8)],[yy(4) yy(8)],[zz(4) zz(8)],'Color','k');
    plot3(xx, yy, zz, 'k.')
    hold on
    
    % FACE 1
    xp  = xx(1:4); yp  = yy(1:4); zp  = zz(1:4); cp  = cc(1:4);
    p   = patch(xp, yp, zp, cp);
    
    % FACE 2
    xp  = xx(5:8); yp  = yy(5:8); zp  = zz(5:8); cp  = cc(5:8);
    p   = patch(xp, yp, zp, cp);
    
    % FACE 3
    xp  = [xx(1) xx(5) xx(8) xx(4)];
    yp  = [yy(1) yy(5) yy(8) yy(4)];
    zp  = [zz(1) zz(5) zz(8) zz(4)];
    cp  = [cc(1) cc(5) cc(8) cc(4)];
    p   = patch(xp, yp, zp, cp);
    
    % FACE 4
    xp  = [xx(2) xx(3) xx(7) xx(6)];
    yp  = [yy(2) yy(3) yy(7) yy(6)];
    zp  = [zz(2) zz(3) zz(7) zz(6)];
    cp  = [cc(2) cc(3) cc(7) cc(6)];
    p   = patch(xp, yp, zp, cp);
    
    % FACE 5
    xp  = [xx(2) xx(6) xx(5) xx(1)];
    yp  = [yy(2) yy(6) yy(5) yy(1)];
    zp  = [zz(2) zz(6) zz(5) zz(1)];
    cp  = [cc(2) cc(6) cc(5) cc(1)];
    p   = patch(xp, yp, zp, cp);
    
    % FACE 6
    xp  = [xx(3) xx(7) xx(8) xx(4)];
    yp  = [yy(3) yy(7) yy(8) yy(4)];
    zp  = [zz(3) zz(7) zz(8) zz(4)];
    cp  = [cc(3) cc(7) cc(8) cc(4)];
    p   = patch(xp, yp, zp, cp);
end
axis equal tight off
colorbar horz
caxis([-400 400])