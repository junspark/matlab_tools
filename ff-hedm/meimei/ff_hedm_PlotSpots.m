clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CRYSTAL SYMMETRY
qsym    = CubSymmetries;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ff_HEDM DATA ACQUISITION PARAMETERS
Dsam        = 910903.504046489;
ome_ini         = 180;
ome_fin         = -180;
ome_step        = -0.5;

num_img_per_layer       = 3;
num_ome_per_img_stack   = 240;

pfimage{1}  = 'O:\meimei_aug14\ge3\hedm_HTUPS_irr_s1_02546.ge3';
pfimage{2}  = 'O:\meimei_aug14\ge3\hedm_HTUPS_irr_s1_02547.ge3';
pfimage{3}  = 'O:\meimei_aug14\ge3\hedm_HTUPS_irr_s1_02548.ge3';
pfdarak     = 'O:\meimei_aug14\ge3\dark_02587.ge3';

pname_ome   = 'O:\li_aug14\hedm_HTUPS_irr_s1_MultiRing\Layer1_ring1_t150_ring2_t50\Output';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ff_HEDM ANALYSIS RESULT
grains	= load('O:\li_aug14\hedm_HTUPS_irr_s1_MultiRing/Layer1_ring1_t150_ring2_t50\Grains.csv');

CENX    = 1021.4810207067;
CENY    = 1039.6284779322;

%%% THRESHOLDING BY COMPLETENESS
Thresh_Completeness = 1.0;
idx_Completeness    = grains(:,24) >= Thresh_Completeness;

%%% THRESHOLDING BY MEAN RADIUS
Thresh_MeanRadius   = 100;
idx_MeanRadius      = grains(:,23) >= Thresh_MeanRadius;

% grains  = grains(idx_MeanRadius, :);
grains  = grains(idx_Completeness, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ROI SIZE (HALF THE WINDOW)
window_size         = 15;
omega_wiggle_step   = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIT OUTPUT
pfname_output   = 'O:\li_aug14\hedm_HTUPS_irr_s1_MultiRing/Layer1_ring1_t150_ring2_t50/Layer1_ring1_t150_ring2_t50_SpotFits.csv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_grains      = size(grains, 1);

num_ome_steps   = (ome_fin - ome_ini)/ome_step;
ome_grid        = linspace(ome_ini + ome_step/2, ome_fin - ome_step/2, num_ome_steps);

num_total_frames    = num_img_per_layer * num_ome_per_img_stack;

frame_number    = 1:1:num_ome_per_img_stack;
frame_number    = repmat(frame_number, 1, num_img_per_layer);

im_stack_number = zeros(1,num_total_frames);
for i = 1:1:num_img_per_layer
    idx1    = (i - 1) * num_ome_per_img_stack + 1;
    idx2    = i * num_ome_per_img_stack;
    im_stack_number(idx1:idx2)  = ones(1,num_ome_per_img_stack) * i;
end

RESRF2APS   = RMatOfQuat(QuatOfESRF2APS);
RMats       = zeros(3, 3, num_grains);
for i = 1:1:num_grains
    RMats(:,:,i)    =  RESRF2APS * reshape(grains(i,2:10), 3, 3)';
end
quat    = ToFundamentalRegionQ(QuatOfRMat(RMats), qsym);
rod     = RodOfQuat(quat);

darkdata    = NreadGE(pfdarak, 1);
data2write  = [];
for i = 1:1:num_grains
    figure(1); clf
    subplot(1,2,1)
    scatter3(grains(:,11), grains(:,12), grains(:,13), 'b', 'o')
    hold on
    scatter3(grains(i,11), grains(i,12), grains(i,13), 'filled', 'r', 'o')
    xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
    axis equal tight; axis([-1000 1000 -1000 1000 -200 200])
    
    subplot(1,2,2)
    PlotFRPerimeter('cubic');
    scatter3(rod(1,:), rod(2,:), rod(3,:), 'b', 'o')
    hold on
    scatter3(rod(1,i), rod(2,i), rod(3,i), 'filled', 'r', 'o')
    axis equal
    
    fname_ome   = sprintf('FitBest_%09d.csv', grains(i,1));
    pfname_ome  = fullfile(pname_ome, fname_ome);
    Spots       = importdata(pfname_ome, ' ', 1);
    
    num_spots   = size(Spots.data, 1);
    SpotComposite           = zeros(window_size * 2 + 1, (window_size * 2 + 1) * num_spots);
    SpotCompositeStackSum   = zeros(window_size * 2 + 1, (window_size * 2 + 1) * num_spots);
    for j = 1:1:num_spots
        detX    = Spots.data(j,8)/200;
        detY    = Spots.data(j,9)/200;
        detX    = -detX + CENX;
        detY    = -detY + (2048 - CENY);
        
        ome     = Spots.data(j,10);
        
        [~, idx] = min(abs(ome_grid - ome));
        
        imdata  = rot90(NreadGE(pfimage{im_stack_number(idx)}, frame_number(idx)) - darkdata, 1);
        
        figure(2); clf
        imagesc(double(imdata), [0 200])
        axis square tight
        hold on
        colorbar vert
        plot(CENX, 2048 - CENY, 'k.')
        plot(detX, detY, 'ro')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GO +/- OMEGA STEPS FROM CENTROID OMEGA POSITIONS
        %%% OMEGA WIGGLE IS NOT IMPLEMENTED HERE YET
        %%% ONLY FITTING FOR OMEGA0 WHERE SPOT WAS FOUND
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idxIni  = idx - omega_wiggle_step;
        idxFin  = idx + omega_wiggle_step;
        
        if idxIni < 1
            idxIni  = 1;
        end
        
        if idxFin > num_ome_steps
            idxFin  = num_ome_steps;
        end
        idxGrid = idxIni:1:idxFin;
        
        SpotStack   = zeros(window_size*2 + 1, window_size*2 + 1, length(idxGrid));
        for k = 1:1:length(idxGrid)
            imdata  = rot90(NreadGE(pfimage{im_stack_number(idxGrid(k))}, frame_number(idxGrid(k))) - darkdata, 1);
            
            ri  = round(detY - window_size);
            rf  = round(detY + window_size);
            
            ci  = round(detX - window_size);
            cf  = round(detX + window_size);
            
            if ri < 0
                ri  = 1;
            end
            if ci < 0
                ci = 1;
            end
            if rf > 2048
                rf = 2048;
            end
            if cf > 2048
                cf = 2048;
            end
            
            xWindow = (ci:cf) - CENX;
            yWindow = (ri:rf) - (2048 - CENY);
            
            [xWindow, yWindow]  = meshgrid(xWindow, yWindow);
            xWindow = xWindow(:);
            yWindow = yWindow(:);
            
            num_window_pixels   = length(yWindow);
            
            rWindow     = sqrt(xWindow(:).*xWindow + yWindow.*yWindow);
            etaWindow   = rad2deg(atan2(yWindow, xWindow)) + 180;
            
            if abs(sum(sign(yWindow))) ~= num_window_pixels && (sum(xWindow < 0) == num_window_pixels)
                idx_eta = yWindow < 0;
                etaShifted  = 360 + etaWindow(idx_eta);
                etaWindow(idx_eta)  = etaShifted;
            end
            
            Ipk = imdata(ri:rf, ci:cf);
            Ipk = double(Ipk(:));
            
            [maxIpk, idxMaxIpk] = max(Ipk);
            p0  = [max(Ipk); 1.75; 0.5; rWindow(idxMaxIpk); etaWindow(idxMaxIpk); 10];
            [p, resnorm]    = lsqcurvefit(@pkGaussian2, p0, [rWindow etaWindow], Ipk);
            
            Ipk_fit0    = pkGaussian2(p0, [rWindow etaWindow]);
            Ipk_fit     = pkGaussian2(p, [rWindow etaWindow]);
            
            tthMean = atand(p(4)*200./Dsam);
            etaMean = p(5);
            
            SpotStack(:,:,k)    = imdata(ri:rf, ci:cf);
            
            figure(2)
            subplot(2,2,1); cla
            imagesc(SpotStack(:,:,k))
            hold on
            axis square tight
            title('measured peak profile ')
            caxis([0 300])
            colorbar vert
            
            subplot(2,2,2); cla
            imagesc(reshape(Ipk_fit0, window_size*2+1, window_size*2+1))
            hold on
            axis square tight
            title('peak profile from initial guess')
            
            subplot(2,2,3); cla
            imagesc(reshape(Ipk_fit, window_size*2+1, window_size*2+1))
            hold on
            axis square tight
            title('peak profile from fit')
            
            subplot(2,2,4); cla
            imagesc(SpotStack(:,:,k) -  reshape(Ipk_fit, window_size*2+1, window_size*2+1))
            hold on
            axis square tight
            title('fit - measured peak')
            caxis([-300 300])
            colorbar vert
            
            ci  = 1 + (window_size * 2 + 1) * (j - 1);
            cf  = (window_size * 2 + 1) * j;
            
            SpotCompositeStackSum(:, ci:cf) = SpotCompositeStackSum(:, ci:cf) + SpotStack(:,:,k);
            if idxGrid(k) == idx
                ci  = 1 + (window_size * 2 + 1) * (j - 1);
                cf  = (window_size * 2 + 1) * j;
                SpotComposite(:, ci:cf) = SpotStack(:,:,k);
                %%% grain# spot# tth eta Amp gamma_r gamma_eta r_fit eta_fit bck grain radius
                data2write  = [data2write; i j tthMean etaMean p' resnorm grains(i,23)];
            end
        end
    end
    return
end

SpotStack
SpotCompositeStackSum
csvwrite(pfname_output, data2write)
