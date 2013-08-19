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
%%% DIC IMAGE SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname   = 'W:\Miller_June13\DIC';
froot   = 'Ti64_State4_';
fext    = 'tif';
ndigits = 5;
fnum    = [1142 1147 1148 1149 1150 1151 1153 1154 1155 1156 1157 1158 1159 1160 1161 1162 1163];
flist   = cell(length(fnum),1);
padding = '0';
pix2mm  = 0.002;    %%% mm / pixel
gauge_length	= 0.354*25.4;
stress  = [0 100 200 300 400 500 600 650 675 700 725 750 775 800 825 850 875];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:length(fnum)
    if fnum(i) < 0
        error('file number must be larger than or equal to 0')
    else
        flist{i,1}  = [froot, ...
            num2str(stress(i)), 'MPa_', ...
            sprintf(['%0', num2str(ndigits), 'd'], fnum(i)), ...
            '.', fext];
    end
end

stress_pydata   = [];
strain_pydata   = [];
for i = 2:1:length(fnum)
    fname_pydata    = [froot, 'stress_strain_', num2str(stress(i-1)), '_', num2str(stress(i)), '.pydata'];
    pfname_pydata   = fullfile('X:\Miller_Jun13', fname_pydata);
    if exist(pfname_pydata, 'file')
        pydata  = load(pfname_pydata);
    end
    stress_pydata   = [stress_pydata; pydata(:,1)];
    strain_pydata   = [strain_pydata; pydata(:,3)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINE ROI & CONTROL POINTS
%%% Generate a regular grid of points, with spacing delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname0  = flist{i,1};   %%% LOAD INITIAL STATE DIC IMAGE
pfname0 = fullfile(pname, fname0);
imdata0 = imread(pfname0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ROI
% Hctr    = 1238;     % Center point in horizontal direction
% Vctr    = 1163;     % Center point in vertical direction
% Hctr    = 1249;     % Center point in horizontal direction
% Vctr    = 908;      % Center point in vertical direction
%%%
Hctr    = 1223;     % Center point in horizontal direction
Vctr    = 1031;      % Center point in vertical direction
delta_h = 2;        % Spacing between control points
delta_v = 2;        % Spacing between control points
ws_h    = 7;       % The control points go from (-ws,-ws) to (ws,ws)
ws_v    = 7;       % The control points go from (-ws,-ws) to (ws,ws)
%%%
% Hctr    = 1223;     % Center point in horizontal direction
% Vctr    = 1031;      % Center point in vertical direction
% delta_h = 5;        % Spacing between control points
% delta_v = 5;        % Spacing between control points
% ws_h    = 15;       % The control points go from (-ws,-ws) to (ws,ws)
% ws_v    = 20;       % The control points go from (-ws,-ws) to (ws,ws)
%%%
% Hctr    = 1238;     % Center point in horizontal direction
% Vctr    = 1036;      % Center point in vertical direction
% delta_h = 14;        % Spacing between control points
% delta_v = 20;        % Spacing between control points
% ws_h    = 15;       % The control points go from (-ws,-ws) to (ws,ws)
% ws_v    = 40;       % The control points go from (-ws,-ws) to (ws,ws)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

figure(1)
im0zoomed   = imdata0(min(v_im):max(v_im), min(h_im):max(h_im));
subplot(1,2,2)
imagesc(im0zoomed)
colormap(gray)
axis equal tight on
hold off
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
    % plot(samY(1,:), -samY(2,:), 'ko')
    legend('ave-u', 'ave-v', 'Location', 'NorthWest')
    % legend('ave-u', 'ave-v', 'samY', 'Location', 'NorthWest')
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
    
    figure(400)
    subplot(1,2,1)
    imagesc(reshape(disp_h', 2*ws_v + 1, 2*ws_h + 1));
    title('displacement - horizontal')
    axis equal tight off
    colorbar vert
    hold off
    
    subplot(1,2,2)
    imagesc(reshape(disp_v', 2*ws_v + 1, 2*ws_h + 1));
    title('displacement - vertical')
    axis equal tight off
    colorbar vert
    hold off
    % return
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
plot(strain_pydata*10^-6, stress_pydata, 'ro-')
plot(dg(:,6), stress, 'bo-')
xlabel('axial strain (-)')
ylabel('axial stress (MPa)')
legend('displacement / length','DIC strain','extensometer strain','Location','NorthWest')
% legend('displacement / length','extensometer strain','Location','NorthWest')
hold off