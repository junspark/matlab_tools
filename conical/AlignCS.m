close all
clear all
clc

%%%
%%% look at images
% pname   = '~/mnt/s1c/hchoo_feb18/dexela/';
% froot   = 'dark_';
% fnumber = 266;
% 
% froot   = sprintf('%s%06d.tif', froot, fnumber);
% pfname  = fullfile(pname, froot);
% bg  = imread(pfname);
% 
% froot   = 'Al_cs_align_';
% fnumber = 267;
% 
% froot   = sprintf('%s%06d.tif', froot, fnumber);
% pfname  = fullfile(pname, froot);
% im  = imread(pfname, 1);
% img = (im-bg);
% % 
% figure(fnumber)
% % imagesc(img,[1700 2100]); colorbar
% imagesc(img,[-10 50]); colorbar
% % imagesc(log(abs(img)),[-10 900]); colorbar
% % imagesc(log(abs(img))); colorbar
% % imagesc(img); colorbar
% axis equal
% title(fnumber)
% return
%%%

%%% analyse samZ scans
% load background
pname   = '~/mnt/s1c/hchoo_feb18/dexela/';
froot   = 'dark_';
fnumber = 266;

froot   = sprintf('%s%06d.tif', froot, fnumber);
pfname  = fullfile(pname, froot);
bg  = imread(pfname);
    
sd  = 566.255;
a0  = 4.0495;
en  = 55.000;
cx  = 1949.719; % in pixels
cy  = 1560.875; % in pixels

hkl2    = [3 4 8 12];

nhkl    = length(hkl2);
nhkl2   = nhkl*2;
tth     = 2 * asin(12.39854/en./(a0./sqrt(hkl2))/2);
radius  = sd * tan(tth) / 0.075; % radius in pixels
% radius0 = round([-fliplr(radius) radius] + 1024);
radius0 = round([-fliplr(radius) radius]);
dt      = 60;
dr      = 20;

% figure(2)
% plot(1:2048, mean(img(cy-dt:cy+dt,:)), radius0, 2000*ones(nhkl2,1),'o')
% return
% radius  = radius0 - 1024;
% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 11)+1.95;
% imno = 940:950;
% % horizontal - beta - rotation about y
% %     0.0077
% % vertical - alpha - rotation about x
% %    -3.6227

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 11)+1.95;
% imno = 951:961;
% % horizontal - beta - rotation about y
% %     2.0722
% % vertical - alpha - rotation about x
% %    -3.5193

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 11)+1.95;
% imno = 973:983;
% % horizontal - beta - rotation about y
% %    -0.3857
% % vertical - alpha - rotation about x
% %    15.1671

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 11)+1.95;
% % imno = 1051:1061;
% imno = 1061:1071;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 51)+1.95;
% imno = 1072:1122;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.915;
% imno = 1123:1148;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.875;
% imno = 1149:1174;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.875;
% imno = 1175:1200;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.875;
% imno = 1201:1226;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.875;
% imno = 1227:1252;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.875;
% imno = 1253:1278;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.875;
% imno = 1279:1304;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 26)+1.875;
% imno = 1305:1330;

% froot = 'testinng_';
% samZ = linspace(-1.0, 1.0, 50)+1.875;
% imno = 1331:1380;

froot = 'testinng_';
samZ = linspace(-0.5, 0.5, 51)+1.875;
imno = 1381:1431;

%%%%%%%%%%%%%%%%%%%%%%%
for mm = length(imno):-1:1
    fnumber = imno(mm);
    fname   = sprintf('%s%06d.tif', froot, fnumber);
    pfname  = fullfile(pname, fname);
    while ~exist(pfname, 'file')
        pause(2)
    end
    
    img = imread(pfname);
    
    figure(1)
    colormap jet
    imagesc(img,[40 2000]); colorbar
    axis equal
    title(imno(mm))

%     pause
    % figure(2)
    % plot(1:3888, mean(img(cy-dt:cy+dt,:)), radius0, ones(nhkl2,1), 'o')
    % plot(radius0, ones(nhkl2,1),'o')
%     pause
    % integrate horizontal sections
    for i=nhkl2:-1:1
%         figure(20)
        radiusv = radius0(i) + cx;
        inth(mm,i) = sum(sum(img(cy-dt:cy+dt, radiusv-dr:radiusv+dr)));
%         imagesc(img(cy-dt:cy+dt, radiusv-dr:radiusv+dr))
%         colorbar
%         title(['hor. pos: ' num2str(radiusv), ' Z = ', num2str(samZ(mm))]);pause(.2)
    end
    
    % integrate vertical sections
    for i=nhkl2:-1:1
        %         figure(20)
        radiush = radius0(i) + cy;
        %         imagesc(img(max([radiush-dr,1]):radiush+dr,cx-dt:cx+dt))
        %         colorbar
        %         title(['vert. pos: ' num2str(radiush), ' Z = ', num2str(samZ(mm))]);pause(.2)
        intv(mm,i) = sum(sum(img(max([radiush-dr,1]):radiush+dr,cx-dt:cx+dt)));
        %         pause(0.5)
    end
end

%%% look at samZ peaks
figure(30)
% contour(log(inth), 40)
contour((inth), 30)
hold on
[px py] = gradient(inth,.25,.25);
quiver(px,py)
title('horizontal symmetry')
xlabel('ring number')
ylabel('samZ number')
hold off

figure(31)
% contour(log(intv), 40)
contour((intv), 30)
hold on
[px py] = gradient(intv,.25,.25);
quiver(px,py)
title('vertical symmetry')

xlabel('ring number')
ylabel('samZ number')
hold off
% return
figure(4)
for i=[nhkl nhkl+1]
    plot(samZ, inth(:,i));hold on
end
hold off; grid on

%% fit samZ center positions
samZ =  samZ(:);
figure(4)
for i = nhkl2:-1:1
    y   = inth(:, i);
    
    p(4)    = min(y);
    p(5)    = (y(1) - y(end))/(samZ(1) - samZ(end));
    p(1)    = max(y)-p(4);
    p(2)    = mean(samZ);
    p(3)    = 0.25;
    
    % lb = [-inf x(1) 0 0.0 0];
    % ub = [inf x(end) inf 1 inf];
    opts = optimset('display','off');
    fit = lsqcurvefit(@gauss, p, samZ, y, [], [], opts);
    plot(samZ, gauss(fit, samZ), samZ, y, 'o')
    title(['horiz ' num2str(i)]); pause(.5)
    xlabel('samZ')
    posh(i) = fit(2);
    areah(i) = fit(1);%*fit(3);
end

figure(5)
for i = nhkl2:-1:1
    y = intv(:, i);
    p(4)    = min(y);
    p(5)    = (y(1) - y(end))/(samZ(1) - samZ(end));
    p(1)    = max(y)-p(4);
    p(2)    = mean(samZ);
    p(3)    = 0.25;
    
    % lb = [-inf x(1) 0 0.0 0];
    % ub = [inf x(end) inf 1 inf];
    opts = optimset('display','off');
    fit = lsqcurvefit(@gauss, p, samZ, y, [], [], opts);
    plot(samZ, gauss(fit, samZ), samZ, y, 'o')
    title(['vert ' num2str(i)]); pause(.5)
    xlabel('samZ')
    posv(i) = fit(2);
    areav(i) = fit(1)*fit(3);
end

% radius  = sd * tan(tth) / 0.2;
% radius  = round([-fliplr(radius) radius] + 1024);
% 
% figure(6)
% plot(radius(nhkl+1:nhkl2)-cx, posh(nhkl+1:nhkl2)-posh(nhkl:-1:1),'-bo'); hold on
% plot(radius(nhkl+1:nhkl2)-cy, posv(nhkl+1:nhkl2)-posv(nhkl:-1:1),'-rx'); hold off
% xlabel('radius [pix]'), ylabel('samZ [mm]'), title('center samZ')
% legend('horizontal', 'vertical')

% hkl2    = [3 4 8 12];
conical_radii = [4.8378, 5.5927, 7.9465, 9.7787];
% conical_radii = [7.9465, 12.4080, 14.03];

disp(asin((posh(nhkl+1:nhkl2)-posh(nhkl:-1:1))./conical_radii))
disp(asin((posv(nhkl+1:nhkl2)-posv(nhkl:-1:1))./conical_radii))

disp(mean(asin((posh(nhkl+1:nhkl2)-posh(nhkl:-1:1))./conical_radii)))
disp(mean(asin((posv(nhkl+1:nhkl2)-posv(nhkl:-1:1))./conical_radii)))

disp('horizontal - beta - rotation about y')
disp(-mean(asin((posh(nhkl+1:nhkl2)-posh(nhkl:-1:1))./conical_radii))*1000*43/1000)
disp('vertical - alpha - rotation about x')
disp(mean(asin((posv(nhkl+1:nhkl2)-posv(nhkl:-1:1))./conical_radii))*1000*148/1000)

% figure(7)
% plot(radius(nhkl+1:nhkl2)-cx, areah(nhkl+1:nhkl2)./areah(nhkl:-1:1),'-bo'); hold on
% plot(radius(nhkl+1:nhkl2)-cy, areav(nhkl+1:nhkl2)./areav(nhkl:-1:1),'-rx'); hold off
% xlabel('radius [pix]'), ylabel('ratio'), title('area ratio')
% legend('horizontal', 'vertical')

[mean(posh) mean(posv) mean([posh posv]) std([posh posv])]

% figure(1000)
% plot(radius-cx, posh,'-bo'); hold on
% plot(radius-cy, posv,'-rx'); hold off
% xlabel('radius [pix]'), ylabel('samZ [mm]'), title('center samZ')
% legend('horizontal', 'vertical')
