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
load('MLS.Stress.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPTIMIZE E / nu
for j = 1:1:nr
    for i = 1:1:na       
        %%%%%%%%%%%%%%%%%%%%%%
        %%% hkl        
        q110    = LS110{i,j}(:,1:3);
        LS110ij = LS110{i,j}(:,4);
        
        sin2psi110 = 1 - q110(:,3).^2;
        
        q112    = LS112{i,j}(:,1:3);
        LS112ij = LS112{i,j}(:,4);
        
        sin2psi112 = 1 - q112(:,3).^2;
        
        idx110 = sin2psi110 < 0.6;
        idx112 = sin2psi112 < 0.6;
        
        p110ft2{i,j}    = polyfit(sin2psi110(idx110), LS110ij(idx110), 1);       
        m110(i,j)       = p110ft2{i,j}(1);
        b110(i,j)       = p110ft2{i,j}(2);
        
        p112ft2{i,j}    = polyfit(sin2psi112(idx112), LS112ij(idx112), 1);       
        m112(i,j)       = p112ft2{i,j}(1);
        b112(i,j)       = p112ft2{i,j}(2);
    end
end

LB      = [25000 0.00];
UB      = [200000 0.50];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 110
%%%%%%%%%%%%%%%%%%%%%%%%%%
S12_110i    = [100000 0.3]';               %%% PROPERTY INITIAL GUESS
%%% ALL DATA
xdata_a = [m110(:) b110(:)];            
ydata_a = [sx(:); sy(:)];
[S12_110f_a,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty, S12_110i, xdata_a, ydata_a, ...
    LB, UB);
yfunc0_a    = OptimizeProperty(S12_110i,  xdata_a);
yfunc_a     = OptimizeProperty(S12_110f_a,xdata_a);

figure,
x   = [r_dvc r_dvc r_dvc r_dvc]';
plot(x, ydata_a, 'k.')
hold on
plot(x, yfunc0_a, 'b.')
plot(x, yfunc_a, 'r.')
title(['\{110\} - A'])

%%% SIGXX - A0
xdata_b = [m110(1,:); b110(1,:)]';
ydata_b = [sx(1,:)]';
[S12_110f_b,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_xx, S12_110i, xdata_b, ydata_b, ...
    LB, UB);
yfunc0_b    = OptimizeProperty_xx(S12_110i,  xdata_b);
yfunc_b     = OptimizeProperty_xx(S12_110f_b,xdata_b);

figure,
x   = r_dvc';
plot(x, ydata_b, 'k.')
hold on
plot(x, yfunc0_b, 'b.')
plot(x, yfunc_b, 'r.')
title(['\{110\} - B'])

%%% SIGXX - A90
xdata_c = [m110(2,:); b110(2,:)]';
ydata_c = [sx(2,:)]';
[S12_110f_c,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_xx, S12_110i, xdata_c, ydata_c, ...
    LB, UB);
yfunc0_c    = OptimizeProperty_xx(S12_110i,  xdata_c);
yfunc_c     = OptimizeProperty_xx(S12_110f_c,xdata_c);

figure,
x   = r_dvc';
plot(x, ydata_c, 'k.')
hold on
plot(x, yfunc0_c, 'b.')
plot(x, yfunc_c, 'r.')
title(['\{110\} - C'])

%%% SIGYY - A0
xdata_d = [m110(1,:); b110(1,:)]';
ydata_d = [sy(1,:)]';
[S12_110f_d,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_yy, S12_110i, xdata_d, ydata_d, ...
    LB, UB);
yfunc0_d    = OptimizeProperty_yy(S12_110i,  xdata_d);
yfunc_d     = OptimizeProperty_yy(S12_110f_d,xdata_d);

figure,
x   = r_dvc';
plot(x, ydata_d, 'k.')
hold on
plot(x, yfunc0_d, 'b.')
plot(x, yfunc_d, 'r.')
title(['\{110\} - D'])

%%% SIGYY - A90
xdata_e = [m110(2,:); b110(2,:)]';
ydata_e = [sy(2,:)]';
[S12_110f_e,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_yy, S12_110i, xdata_e, ydata_e, ...
    LB, UB);
yfunc0_e    = OptimizeProperty_yy(S12_110i,  xdata_e);
yfunc_e     = OptimizeProperty_yy(S12_110f_d,xdata_e);

figure,
x   = r_dvc';
plot(x, ydata_e, 'k.')
hold on
plot(x, yfunc0_e, 'b.')
plot(x, yfunc_e, 'r.')
title(['\{110\} - E'])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 112
%%%%%%%%%%%%%%%%%%%%%%%%%%
S12_112i    = [100000 0.3]';               %%% PROPERTY INITIAL GUESS
%%% ALL DATA
xdata_a = [m112(:) b112(:)];
ydata_a = [sx(:); sy(:)];
[S12_112f_a,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty, S12_112i, xdata_a, ydata_a, ...
    LB, UB);
yfunc0_a    = OptimizeProperty(S12_112i,  xdata_a);
yfunc_a     = OptimizeProperty(S12_112f_a,xdata_a);
disp(S12_112f_a)
figure,
plot(ydata_a, 'k.')
hold on
plot(yfunc0_a, 'b.')
plot(yfunc_a, 'r.')
title(['\{112\} - A'])

%%% SIGXX - A0
xdata_b = [m112(1,:); b112(1,:)]';      
ydata_b = [sx(1,:)]';
[S12_112f_b,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_xx, S12_112i, xdata_b, ydata_b, ...
    LB, UB);
yfunc0_b    = OptimizeProperty_xx(S12_112i,  xdata_b);
yfunc_b     = OptimizeProperty_xx(S12_112f_b,xdata_b);
disp(S12_112f_b)
figure,
plot(ydata_b, 'k.')
hold on
plot(yfunc0_b, 'b.')
plot(yfunc_b, 'r.')
title(['\{112\} - B'])

%%% SIGXX - A90
xdata_c = [m112(2,:); b112(2,:)]';
ydata_c = [sx(2,:)]';
[S12_112f_c,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_xx, S12_112i, xdata_c, ydata_c, ...
    LB, UB);
yfunc0_c    = OptimizeProperty_xx(S12_112i,  xdata_c);
yfunc_c     = OptimizeProperty_xx(S12_112f_c,xdata_c);
disp(S12_112f_c)
figure,
plot(ydata_c, 'k.')
hold on
plot(yfunc0_c, 'b.')
plot(yfunc_c, 'r.')
title(['\{112\} - C'])

%%% SIGYY - A0
xdata_d = [m112(1,:); b112(1,:)]';
ydata_d = [sy(1,:)]';
[S12_112f_d,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_yy, S12_112i, xdata_d, ydata_d, ...
    LB, UB);
yfunc0_d    = OptimizeProperty_yy(S12_112i,  xdata_d);
yfunc_d     = OptimizeProperty_yy(S12_112f_d,xdata_d);
disp(S12_112f_d)
figure,
plot(ydata_d, 'k.')
hold on
plot(yfunc0_d, 'b.')
plot(yfunc_d, 'r.')
title(['\{112\} - D'])

%%% SIGYY - A90
xdata_e = [m112(2,:); b112(2,:)]';      
ydata_e = [sy(2,:)]';
[S12_112f_e,~,~,ef,output,lambda,~] = lsqcurvefit(@OptimizeProperty_yy, S12_112i, xdata_e, ydata_e, ...
    LB, UB);
yfunc0_e    = OptimizeProperty_yy(S12_112i,  xdata_e);
yfunc_e     = OptimizeProperty_yy(S12_112f_d,xdata_e);
disp(S12_112f_e)
figure,
plot(ydata_e, 'k.')
hold on
plot(yfunc0_e, 'b.')
plot(yfunc_e, 'r.')
title(['\{112\} - E'])

disp([ ...
    round(S12_110f_a(1)/1000) S12_110f_a(2);
    round(S12_110f_b(1)/1000) S12_110f_b(2);
    round(S12_110f_c(1)/1000) S12_110f_c(2);
    round(S12_110f_d(1)/1000) S12_110f_d(2);
    round(S12_110f_e(1)/1000) S12_110f_e(2);
    round(S12_112f_a(1)/1000) S12_112f_a(2);
    round(S12_112f_b(1)/1000) S12_112f_b(2);
    round(S12_112f_c(1)/1000) S12_112f_c(2);
    round(S12_112f_d(1)/1000) S12_112f_d(2);
    round(S12_112f_e(1)/1000) S12_112f_e(2);
    ])