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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIT LINE THROUGH SIN2PSI VS. STRAIN
x   = 0:0.1:1;
for j = 1:1:nr
    fig = figure(j);
    orient rotated
    set(fig, 'position', [1681 27 1280 940], ...
        'PaperOrientation', 'portrait')
    for i = 1:1:na
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%% 110
        q   = LS110{i,j}(:,1:3);
        LS  = LS110{i,j}(:,4);
        sin2psi = 1 - q(:,3).^2;
        
        p110{i,j}   = polyfit(sin2psi, LS, 1);
        LSfit       = polyval(p110{i,j}, x);
        p110ft{i,j} = polyfit(sin2psi(1:4:end), LS(1:4:end), 1);
        LSfitft     = polyval(p110ft{i,j}, x);
        
        idx = sin2psi < 0.6 & sin2psi > 0.05;
        p110ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);
        LSfitft2        = polyval(p110ft2{i,j}, x);
        
        figure(j)
        set(subplot(2,2,1+(i-1)*2), ...
            'FontSize', 14, ...
            'FontWeight', 'bold')
        hold on
        plot(sin2psi(idx), LS(idx), 'ko', ...
            'MarkerSize', 6, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'k')
%         plot(x, LSfit, 'b-', ...
%             'LineWidth', 3)
%         plot(x, LSfitft, 'r--', ...
%             'LineWidth', 3)
        plot(x, LSfitft2, 'k-', ...
            'LineWidth', 3)
        axis([0 1 LSRange])
        title(['\{11.0\}'])
        %xlabel('\sin^2 \Psi', 'interpreter', 'tex')
        ylabel('lattice strain')
        grid on
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%% 112
        q   = LS112{i,j}(:,1:3);
        LS  = LS112{i,j}(:,4);
        sin2psi = 1 - q(:,3).^2;
        
        p112{i,j}   = polyfit(sin2psi, LS, 1);
        LSfit       = polyval(p112{i,j}, x);
        p112ft{i,j} = polyfit(sin2psi(1:4:end), LS(1:4:end), 1);
        LSfitft     = polyval(p112ft{i,j}, x);
        
        idx = sin2psi < 0.6 & sin2psi > 0.05;
        p112ft2{i,j}    = polyfit(sin2psi(idx), LS(idx), 1);
        LSfitft2        = polyval(p112ft2{i,j}, x);
        
        figure(j)
        set(subplot(2,2,2+(i-1)*2), ...
            'FontSize', 14, ...
            'FontWeight', 'bold')
        hold on
        plot(sin2psi(idx), LS(idx), 'ko', ...
            'MarkerSize', 6, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'k')
%         plot(x, LSfit, 'b-', ...
%             'LineWidth', 3)
%         plot(x, LSfitft, 'r--', ...
%             'LineWidth', 3)
        plot(x, LSfitft2, 'k-', ...
            'LineWidth', 3)
        axis([0 1 LSRange])
        title(['\{11.2\}'])
        %xlabel('\sin^2 \Psi', 'interpreter', 'tex')
        ylabel('lattice strain')
        grid on
    end
    fig_fname   = ['sin2psi.fit.pos.', num2str(j), '.png'];
    disp(fig_fname)
    %return
    saveas(fig, fig_fname, 'png')
    pause(1)
    close(fig)
end