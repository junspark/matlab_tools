clear all
close all
clc

addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CRYSTAL SYMMETRY
wsname  = 'wscub3x';      % workspace name

pfname_midas_ps = '/home/beams/S1IDUSER/mnt/orthros/faber_mar21_midas/ff/s9/faber_sam1_s9_crack_ff2_Layer81_Analysis_Time_2021_03_07_07_25_19/ps_sam1_s9.txt';

%%% NOTE LAYER NUMBER AND pfimage NEED TO MATCH
pname_grains	= '/home/beams/S1IDUSER/mnt/orthros/faber_mar21_midas/ff/s9/faber_sam1_s9_crack_ff2_Layer81_Analysis_Time_2021_03_07_07_25_19';
pfimage         = '/home/beams/S1IDUSER/mnt/orthros/faber_mar21_data/faber_mar21/ge3/faber_sam1_s9_crack_ff2_000133.ge3';
pfdark          = '/home/beams/S1IDUSER/mnt/orthros/faber_mar21_data/faber_mar21/ge3/dark_before_000052.ge3';

%%% FIT OUTPUT
pfname_output   = './hedm_HTUPS_irr_s1_MultiRing/Layer1_ring1_t150_ring2_t50_SpotFits.csv';

%%%% ff_HEDM ANALYSIS RESULT
pname_ome       = 'O:\faber_mar21_midas\ff\s9\faber_sam1_s9_crack_ff2_Layer1_Analysis_Time_2021_03_07_00_17_55\Output';

%%% THRESHOLDING BY COMPLETENESS; IGNORE GRAINS WITH LOWER COMPLETENESS
Thresh_Completeness = 0.7;

%%% THRESHOLDING BY MEAN RADIUS; IGNORE GRAINS WITH SMALLER RADIUS
Thresh_MeanRadius   = 100;

%%% ROI SIZE (HALF THE WINDOW)
window_size         = 3;    %%% PIXEL
omega_wiggle_step   = 2;    %%% FRAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load workspace for fundamental region.
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname)

% LOAD MIDAS PARAMETERS
midas_params = parseMIDASparams(pfname_midas_ps, pfimage);

% LOAD GRAIN INFORMATION
%%% FOR THE FULL STACK - IGNORE
% for iii = 1:1:length(pname_grains)
%     pfname_grains{iii}  = fullfile(pname_grains{iii}, 'Grains.csv');
% end
% Grains      = parseGrainData(pfname_grains, ws.frmesh.symmetries, ...
%     'CrdSystem', 'APS', ...
%     'LabToSample', 0, ...
%     'C_xstal', BuildElasticityMatrix([386 120 132]*1000, 'Symmetry', 'cubic'));
Grains_midas    = parseGrainData_OneLayer(pname_grains, ws.frmesh.symmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', BuildElasticityMatrix([386 120 132]*1000, 'Symmetry', 'cubic'), ...
    'OutputReflectionTable', true);
Grains_midas.analysis_params    = midas_params;

% SOME FILTERING
idx_Completeness    = [Grains_midas.grains.Completeness] >= Thresh_Completeness;
idx_MeanRadius      = [Grains_midas.grains.GrainRadius] >= Thresh_MeanRadius;
idx_filtered        = idx_Completeness & idx_MeanRadius;

% LOAD DARK
numframes_dark  = CalcNumFramesGE(pfdark);
for iii = 1:1:numframes_dark
    if iii == 1
        bkgpattern  = NreadGE(pfdark, iii);
    else
        bkgpattern  = bkgpattern + NreadGE(pfdark, iii);
    end
end
bkgpattern  = bkgpattern./numframes_dark;

data2write  = [];

com = [Grains_midas.grains.COM];
rod = [Grains_midas.grains.rod];
figure(1); clf
subplot(1,2,1)
scatter3(com(1,:), com(2,:), com(3,:), 'b', 'o')
xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
axis equal tight; axis([-1000 1000 -1000 1000 -200 200])
hold on

subplot(1,2,2)
PlotFRPerimeter('cubic');
scatter3(rod(1,:), rod(2,:), rod(3,:), 'b', 'o')
hold on; axis equal

for iii = 1:1:Grains_midas.nGrains
    if idx_filtered(iii)
        disp(sprintf('grain %d - meets filters; watching ...', Grains_midas.grains(iii).GrainID))
        
        figure(1); subplot(1,2,1)
        scatter3(com(1,iii), com(2,iii), com(3,iii), 'filled', 'r', 'o')
        
        figure(1); subplot(1,2,2)
        scatter3(rod(1,iii), rod(2,iii), rod(3,iii), 'filled', 'r', 'o')
        
        %%% MAKE SPOT COMPOSITE & EXTRACT SPOT FIT DATA
        num_spots   = Grains_midas.grains(iii).ReflectionTable.num_spots;
        SpotComposite   = zeros(window_size * 2 + 1, (window_size * 2 + 1) * num_spots);
        
        ReflectionTable     = Grains_midas.grains(iii).ReflectionTable;
        Spots_csv_DetHCrd   = ReflectionTable.Spots_csv_DetHCrd;
        Spots_csv_DetVCrd   = ReflectionTable.Spots_csv_DetVCrd;
        
        figure(2); clf
        scatter(Spots_csv_DetHCrd, Spots_csv_DetVCrd, 'k', 'o')
        hold on; axis equal; grid on
        axis([0 2048 0 2048])
        for jjj = 1:1:num_spots
            figure(2)
            scatter(Spots_csv_DetHCrd(jjj), Spots_csv_DetVCrd(jjj), 'filled', 'k', 'o')
            
            ome_average             = (ReflectionTable.RingNr_csv_ome_min(jjj) + ReflectionTable.RingNr_csv_ome_max(jjj))/2;
            [~, spot_frame_number]  = min(abs(midas_params.ome_grid - ome_average));
            
            %%% PATTERN AT THE CENTER OF OMEGA RANGE
            pattern_data        = rot90(NreadGE(pfimage, spot_frame_number + midas_params.junk_frames) - bkgpattern, 1);
            
            figure(3); clf
            imagesc(pattern_data)
            colorbar vert; colormap jet; caxis([0 50])
            axis equal tight; hold on; grid on
            plot(midas_params.CENX, 2048 - midas_params.CENY, 'kx')  %%% OK
            plot(Spots_csv_DetHCrd(jjj), 2048 - Spots_csv_DetVCrd(jjj), 'ko')
            
            % pattern_sum_data    = rot90(ReadSUM('/home/beams/S1IDUSER/mnt/orthros/faber_mar21_bc/faber_sam1_s9_crack_ff2/ge3/faber_sam1_s9_crack_ff2_000133.ge3.sum'));
            
            % figure(4); clf
            % imagesc(pattern_sum_data)
            % colorbar vert; colormap jet; caxis([-500 10000])
            % axis equal tight; hold on; grid on
            % plot(midas_params.CENX, 2048 - midas_params.CENY, 'kx')
            % plot(Spots_csv_DetHCrd(jjj), 2048 - Spots_csv_DetVCrd(jjj), 'ko')
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GO +/- OMEGA STEPS FROM CENTROID OMEGA POSITIONS
            %%% OMEGA WIGGLE IS NOT IMPLEMENTED HERE YET
            %%% ONLY FITTING FOR OMEGA0 WHERE SPOT WAS FOUND
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            idxIni  = spot_frame_number - omega_wiggle_step;
            idxFin  = spot_frame_number + omega_wiggle_step;
            
            if idxIni < 1
                idxIni  = 1;
            end
            
            if idxFin > midas_params.num_frames
                idxFin  = midas_params.num_frames;
            end
            idxGrid = idxIni:1:idxFin;
            
            SpotStack   = zeros(window_size*2 + 1, window_size*2 + 1, length(idxGrid));
            detX    = Spots_csv_DetHCrd(jjj);
            detY    = 2048 - Spots_csv_DetVCrd(jjj);
            for kkk = 1:1:length(idxGrid)
                pattern_data    = rot90(NreadGE(pfimage, idxGrid(kkk) + midas_params.junk_frames) - bkgpattern, 1);
                
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
                
                SpotStack(:,:,kkk)    = pattern_data(ri:rf, ci:cf);
            end
            
            for kkk = 1:1:length(idxGrid)
                figure(100+kkk); clf
                imagesc(SpotStack(:,:,kkk))
                colorbar vert; colormap jet; caxis([0 max(SpotStack(:))])
                axis equal tight; grid on
                title(sprintf('measured spot at omega = %2.1f deg', midas_params.ome_grid(idxGrid(kkk))))
            end
            pause
        end
        return
    else
        disp(sprintf('grain %d - does not meet filters; skipping ...', Grains_midas.grains(iii).GrainID))
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return



for i = 1:1:num_grains
    
    
    
    
    
    
    
    
    
    
    for j = 1:1:num_spots
        detX    = Spots.data(j,8)/200;
        detY    = Spots.data(j,9)/200;
        detX    = -detX + CENX;
        detY    = -detY + (2048 - CENY);
        
        ome     = Spots.data(j,10);
        
        [~, idx] = min(abs(ome_grid - ome));
        
        pattern_data  = rot90(NreadGE(pfimage{im_stack_number(idx)}, frame_number(idx)) - bkgpattern, 1);
        
        figure(2); clf
        imagesc(double(pattern_data), [0 200])
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
            pattern_data  = rot90(NreadGE(pfimage{im_stack_number(idxGrid(k))}, frame_number(idxGrid(k))) - bkgpattern, 1);
            
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
            
            Ipk = pattern_data(ri:rf, ci:cf);
            Ipk = double(Ipk(:));
            
            [maxIpk, idxMaxIpk] = max(Ipk);
            p0  = [max(Ipk); 1.75; 0.5; rWindow(idxMaxIpk); etaWindow(idxMaxIpk); 10];
            [p, resnorm]    = lsqcurvefit(@pkGaussian2, p0, [rWindow etaWindow], Ipk);
            
            Ipk_fit0    = pkGaussian2(p0, [rWindow etaWindow]);
            Ipk_fit     = pkGaussian2(p, [rWindow etaWindow]);
            
            tthMean = atand(p(4)*200./Dsam);
            etaMean = p(5);
            
            SpotStack(:,:,k)    = pattern_data(ri:rf, ci:cf);
            
            figure(2)
            subplot(2,2,1); cla
            imagesc(SpotStack(:,:,k))
            hold on
            axis square tight
            title('measured peak profile ')
            
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
            
            if idxGrid(k) == idx
                ci  = 1 + (window_size * 2 + 1) * (j - 1);
                cf  = (window_size * 2 + 1) * j;
                SpotComposite(:, ci:cf) = SpotStack(:,:,k);
                %%% grain# spot# tth eta Amp gamma_r gamma_eta r_fit eta_fit bck grain radius
                data2write  = [data2write; i j tthMean etaMean p' resnorm grains(i,23)];
            end
        end
    end
end

csvwrite(pfname_output, data2write)
