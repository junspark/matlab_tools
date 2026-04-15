clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES
% 1. ASSUMES THAT THE Y DIRECTION (APS CRD) IS HORIZONTAL DIRECTION IN THE DIC
% IMAGE
% 2. EXAMINE CAREFULLY THE COORDINATE SYSTEM(S) BEFORE USING THIS SCRIPT WITH
% IN-SITU LOADING
% 3. EXAMPLE IMAGES PROVIDED BY PROF. A. BEAUDOIN AT UIUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXPERIMENT SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_info_file = 'dic_exp_info.md';     % path to experiment config file

[pname, froot, fext, ndigits, fini, ffin, finc, ...
    delta_h, delta_v, Hctr, Vctr, ws_h, ws_v, pix2mm, gauge_length] = dic_read_exp_setup(exp_info_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILD FILE LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnum  = fini:finc:ffin;
flist = cell(length(fnum), 1);

for i = 1:1:length(fnum)
    if fnum(i) < 0
        error('file number must be larger than or equal to 0')
    else
        flist{i,1} = sprintf('%s_%s.%s', froot, ...
            sprintf(['%0', num2str(ndigits), 'd'], fnum(i)), ...
            fext);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINE ROI & CONTROL POINTS
%%% Generate a regular grid of points, with spacing delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname0  = flist{1,1};   %%% LOAD INITIAL STATE DIC IMAGE
pfname0 = fullfile(pname, fname0);
imdata0 = imread(pfname0);

h_im    = -(ws_h*delta_h):delta_h:(ws_h*delta_h);
v_im    = -(ws_v*delta_v):delta_v:(ws_v*delta_v);
[h_im,v_im] = meshgrid(h_im,v_im);
h_im	= h_im(:) + Hctr;
v_im    = v_im(:) + Vctr;
pts     = [h_im, v_im];

[num_v, num_h]  = size(imdata0);

figure(1)
subplot(1,2,1)
imagesc(imdata0)
colormap(gray)
axis equal tight on
hold on
plot(Hctr, Vctr, 'ro')
plot(h_im, v_im, 'b.')
hold off
title('Full Image')
xlabel('horizontal pixel index')
ylabel('vertical pixel index')

figure(1)
im0zoomed   = imdata0(min(v_im):max(v_im), min(h_im):max(h_im));
subplot(1,2,2)
imagesc(im0zoomed)
colormap(gray)
axis equal tight on
hold off
title('Zoomed ROI')
xlabel('horizontal pixel index')
ylabel('vertical pixel index')
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START DIC OPERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dg      = zeros(length(fnum), 8);               % COEFF ARRAY
fname0  = flist{1,1};                           % LOAD INITIAL STATE DIC IMAGE
pfname0 = fullfile(pname, fname0);

imdata_curr     = imread(pfname0);
imdata_curr_pts = pts;

for i = 2:1:length(fnum)
    imdata_prev     = imdata_curr;              % last image becomes the reference
    imdata_prev_pts = imdata_curr_pts;

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('analyzing %s\n', flist{i,1})
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    % Read in the next image
    pfname      = fullfile(pname, flist{i,1});  %%% LOAD INITIAL STATE DIC IMAGE

    imdata_curr = imread(pfname);

    % Use the cpcorr function of MATLAB to compute the displacement of the
    % points in the image with the cut, relative to base image
    imdata_curr_pts = cpcorr_APS(imdata_curr_pts, imdata_prev_pts, imdata_curr, imdata_prev);

    %%%%%%%%%%%%%%%%%
    % Compute the displacement
    % Begin by setting up the displacement
    % For pixel indices, the row increases downward, while the column
    % increases to the right. Pixel indices are integer values,
    % and range from 1 to the length of the row or column.
    % So, the y-coordinate is flipped.
    disp_h  = (imdata_curr_pts(:,1) - imdata_prev_pts(:,1));
    disp_v  = (-(imdata_curr_pts(:,2) - imdata_prev_pts(:,2)));

    hc  = imdata_prev_pts(:,1);
    vc  = num_v-imdata_prev_pts(:,2);

    %%%%%%%%%%%%%%%%%
    % Develop the displacement field using least squares
    % EQ 9 in Meas. Sci. Technol. 20 (2009) 062001
    A       = [ones(length(hc), 1), hc, vc];
    uCoef   = lsqr(A,disp_h);
    vCoef   = lsqr(A,disp_v);

    % Save the displacement gradient terms, recovered from the least squares fit
    if i == 1
        dg(i,1) = mean(A*uCoef)*pix2mm;
        dg(i,2) = mean(A*vCoef)*pix2mm;

        dg(i,3) = uCoef(2);
        dg(i,4) = uCoef(3);
        dg(i,5) = vCoef(2);
        dg(i,6) = vCoef(3);
    else
        dg(i,1) = dg(i-1,1) + mean(A*uCoef)*pix2mm;
        dg(i,2) = dg(i-1,2) + mean(A*vCoef)*pix2mm;

        dg(i,3) = dg(i-1,3) + uCoef(2);
        dg(i,4) = dg(i-1,4) + uCoef(3);
        dg(i,5) = dg(i-1,5) + vCoef(2);
        dg(i,6) = dg(i-1,6) + vCoef(3);
    end
    %%% DIFFERENCE
    udiff   = mean((A*uCoef - disp_h));
    vdiff   = mean((A*vCoef - disp_v));

    dg(i,7) = udiff;
    dg(i,8) = vdiff;

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Plot displacement with the specimen image on the background
    figure(100)
    imdata_curr_zoomed  = imdata_curr(min(v_im):max(v_im), min(h_im):max(h_im));
    h1  = imagesc(imdata_curr_zoomed);
    axis equal tight on
    colormap(gray)
    hold on;

    ax1 = gca;
    ax2 = axes('position',get(ax1,'position'));
    h2 = quiver(hc, vc, disp_h, disp_v);
    set(ax2,'color','none');
    axis equal tight off
    title('Displacement')
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Plot du/dh, dv/dv, SHEAR, SPIN
    figure(200)
    subplot(1,3,1)
    plot(fnum(1:i), dg(1:i,1), 'bo-')
    hold on
    plot(fnum(1:i), dg(1:i,2), 'go-')
    legend('ave-u', 'ave-v', 'Location', 'NorthWest')
    xlabel('image sequence number')
    ylabel('displacement (mm)')
    hold off

    subplot(1,3,2)
    plot(fnum(1:i), dg(1:i,3), 'bo-')
    hold on
    plot(fnum(1:i), dg(1:i,6), 'go-')
    plot(fnum(1:i), dg(1:i,4) + dg(1:i,5), 'ro-')
    plot(fnum(1:i), dg(1:i,4) - dg(1:i,5), 'co-')
    legend('du/dx','dv/dy','shear','spin','Location','NorthWest')
    xlabel('image sequence number')
    ylabel('partials (-)')
    hold off

    subplot(1,3,3)
    plot(fnum(1:i), dg(1:i,7)*100, 'b-')
    hold on
    plot(fnum(1:i), dg(1:i,8)*100, 'g-')
    legend('u-error','v-error','Location','NorthWest')
    xlabel('image sequence number')
    ylabel('mean %-diff')
    hold off
end

% Plot du/dh, dv/dv, SHEAR, SPIN
figure(200)
subplot(1,3,1)
plot(fnum, dg(:,1), 'bo-')
hold on
plot(fnum, dg(:,2), 'go-')
legend('ave-u', 'ave-v', 'Location', 'NorthWest')
xlabel('image sequence number')
ylabel('displacement (mm)')
hold off

subplot(1,3,2)
plot(fnum, dg(:,3), 'bo-')
hold on
plot(fnum, dg(:,6), 'go-')
plot(fnum, dg(:,4) + dg(1:i,5), 'ro-')
plot(fnum, dg(:,4) - dg(1:i,5), 'co-')
legend('du/dx','dv/dy','shear','spin','Location','NorthWest')
xlabel('image sequence number')
ylabel('partials (-)')
hold off

subplot(1,3,3)
plot(fnum, dg(:,7)*100, 'b-')
hold on
plot(fnum, dg(:,8)*100, 'g-')
legend('u-error','v-error','Location','NorthWest')
xlabel('image sequence number')
ylabel('mean %-diff')

figure(300)
plot(dg(:,2)./gauge_length, stress, 'ko-')
hold on
plot(dg(:,6), stress, 'bo-')
xlabel('axial strain (-)')
ylabel('axial stress (MPa)')
hold off
