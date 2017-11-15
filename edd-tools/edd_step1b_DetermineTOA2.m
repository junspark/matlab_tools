%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Calibrate TOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars all;
close all;
%%% Abraham_oct17
cegv = readedd3_6bm('/home/beams29/S6BM/2017-3/abraham_oct17/cal1_1000x200'); cegv = cegv(1);  % (hv: 1000x200)


load('data/Ce.emission.data');
Bi_emission = [77.1079   74.8148   87.343    10.8388   10.73091  13.0235   12.9799   15.2477];

hc = 12.398419057638671;

Ch         = [1:8192];
ChToEnergy = [0.03477932 0.11038588];     % by Co57 (3 gamma emission), 201710. det-1, 50x200  [0 0.09]
%ChToEnergy = [0.03477582 0.02639247];     % by Co57 (3 gamma emission), 201710. det-2, 50x200  [0 0.025]


E      = Ch*ChToEnergy(1)+ChToEnergy(2);

%%% Oct-17, 1000x200 (10/20/2017)
TOA    = 3;    %(rnx1E3=0.86534750);    %det-1

detno = 1;

% ceria lattice parameters
lat = 5.411651;
[dhkl,hkls] = d0('fcc',lat,15,0);
Ehkl = hc/2/sind(TOA/2)./dhkl;

% plot figure type  [Ch E d S^2 position-intensity]
if_plot = [1 1 0 0 1];   % [raw E d ? z-scan]

% average over position
pos = 1:23;    % Oct.17, det-1 (50x300)

ch_range = [1500:4000];   % CeO2, TOA=3.0, oct.17

cc = colormap(lines(25));

if if_plot(1)
    x = Ch;
    figure(1); clf('reset')
    %set(gcf,'Position',[200 400 800 600])
    line(x,sum(cegv.data{detno}(pos,:),1)/length(pos)/cegv.exp_time,'color',cc(3,:),'displayname',sprintf('CeO2-GV-scan-det%d',detno));
    set(gca,'yscale','log','box','on');
    grid('on');
    xlim([0 4000])
    xlabel('Channel');
    ylabel('Intensity (cps)');
    legend('toggle')
end

if if_plot(2)
    x = E;
    figure(2); clf('reset')
    %set(gcf,'Position',[200 400 800 600])
    line(x,sum(cegv.data{detno}(pos,:),1)/length(pos)/60,'color',cc(3,:),'displayname',sprintf('CeO2-GV-scan-det%d',detno));
    line(Ehkl,ones(length(Ehkl),1),'marker','^','color','b','linestyle','none','displayname','CeO2\_hkl');
    line(Ce_emission(1,:)/1000,Ce_emission(2,:),'marker','^','color','g','linestyle','none','displayname','Ce-emission');
    
    set(gca,'yscale','log','box','on');
    grid('on')
    xlim([20 150])
    xlabel('Energy (KeV)');
    ylabel('Intensity (cps)');
    legend('toggle')
end

if if_plot(5)
    figure(5); clf('reset')
    %set(gcf,'Position',[200 400 800 600])
    set(gca,'box','on','linewidth',2,'fontsize',13)
    x = cegv.motorpos;
    line(x,sum(cegv.data{detno}(:,ch_range),2),'marker','o','color',cc(2,:),'displayname',sprintf('Det-%d',detno))
    %set(gca,'yscale','log');
    xlabel('Position (mm)','fontsize',13);
    ylabel('Intensity (cps)','fontsize',13);
    title(sprintf('Integrated Intensity Ch: %d - %d',ch_range(1),ch_range(end)),'fontsize',17);
    leg = legend('toggle'); set(leg,'fontsize',15);
end

x = E;

% %%
clear p 

% initial guess of peak position (TOA = 3)
pcen0 = {[986 997],[1127 1156], 2174, 2510, 3553, [4167 4356]}; % Oct.17, 1000x200 (det-1)

pk_id = {[1 2] [3 4] 5 6 7 [8 9]};                    % TOA=3

vis_fit = 1;
pos = 1;
sum_range = 1:23;    % Oct.17, Abraham(det-1)

p.cen = zeros(length(pos),length(pk_id));

for j = pos

for i = 1:length(pcen0)
    xdata = min(pcen0{i})-100:max(pcen0{i})+100;

    ydata = sum(cegv.data{detno}(sum_range,xdata),1);   % cegv, sum

    npeaks = length(pcen0{i});
    switch npeaks
        case 1
            fn = {'psv1' 'backg1'}; npar = [4 1];
            pfix = [nan nan nan nan nan]; % fix fitting parameters (currently none are fixed)
            opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
            p02 =  [max(ydata)-min(ydata) pcen0{i}     15 0.5  min(ydata)];
            pUB2 = [       max(ydata)*10  pcen0{i}+10 inf 1.0 +inf];
            pLB2 = [                   0  pcen0{i}-10   0   0 -inf];
            y0   = sumfun1(p02,pfix,npar,fn,xdata);
            
            [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
            yfit = sumfun1(pfit,pfix,npar,fn,xdata);
            txt_cen  = sprintf('Center: %2.1f',pfit(2));
            txt_fwhm = sprintf('FWHM: %2.1f',pfit(3));
            %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(5),pfit(6));
            txt_bg   = sprintf('BG: %2.6f',pfit(5));
            p.cen(j,pk_id{i}) = pfit(2);
        case 2
            fn = {'psv1' 'psv1' 'backg1'}; npar = [4 4 1];
            pfix = [nan nan nan nan nan nan nan nan nan]; % fix fitting parameters (currently none are fixed)
            opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
            p02 =  [max(ydata)-min(ydata) pcen0{i}(1)    15 0.5  max(ydata)/3 pcen0{i}(2)    15 0.5 min(ydata)];
            pUB2 = [       max(ydata)*10  pcen0{i}(1)+5 inf 1.0 max(ydata)*10 pcen0{i}(2)+5 inf 1.0       +inf];
            pLB2 = [                   0  pcen0{i}(1)-5   0   0             0 pcen0{i}(2)-5   0   0       -inf];
            y0   = sumfun1(p02,pfix,npar,fn,xdata);
            
            [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
            yfit = sumfun1(pfit,pfix,npar,fn,xdata);
            txt_cen  = sprintf('Center: %2.1f / %2.1f',pfit(2),pfit(6));
            txt_fwhm = sprintf('FWHM: %2.1f / %2.1f',pfit(3),pfit(7));
            %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(9),pfit(10));
            txt_bg   = sprintf('BG: %2.6f',pfit(9));
            p.cen(j,pk_id{i}) = [pfit(2) pfit(6)];
        case 3
            fn = {'psv1' 'psv1' 'psv1' 'backg1'}; npar = [4 4 4 1];
            pfix = [nan nan nan nan nan nan nan nan nan nan nan nan nan]; % fix fitting parameters (currently none are fixed)
            opt2 = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-6,'tolf',1e-7,'MaxIter',500);
            p02 =  [max(ydata)-min(ydata) pcen0{i}(1)    15 0.5 max(ydata)/3 pcen0{i}(2)    15 0.5  max(ydata)-min(ydata) pcen0{i}(3)    15 0.5 min(ydata)];
            pUB2 = [       max(ydata)*10  pcen0{i}(1)+5 inf 1.0 max(ydata)*10  pcen0{i}(2)+5 inf 1.0 max(ydata)*10 pcen0{i}(3)+5 inf 1.0       +inf];
            pLB2 = [                   0  pcen0{i}(1)-5   0   0             0  pcen0{i}(2)-5   0   0             0 pcen0{i}(3)-5   0   0       -inf];
            y0   = sumfun1(p02,pfix,npar,fn,xdata);
            
            [pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p02,pLB2,pUB2,opt2,pfix,npar,fn,xdata,ydata,ones(size(ydata)));
            yfit = sumfun1(pfit,pfix,npar,fn,xdata);
            txt_cen  = sprintf('Center: %2.1f / %2.1f / %2.1f',pfit(2),pfit(6),pfit(10));
            txt_fwhm = sprintf('FWHM: %2.1f / %2.1f / %2.1f',pfit(3),pfit(7),pfit(11));
            %txt_bg   = sprintf('BG: %2.6f / %2.6f',pfit(9),pfit(10));
            txt_bg   = sprintf('BG: %2.6f',pfit(13));
            p.cen(j,pk_id{i}) = [pfit(2) pfit(6) pfit(10)];
        otherwise
            fprintf('Multi-peaks (>3) fitting not supported yet!!\n');
    end
    
    if vis_fit
    figure(6);
    subplot(4,1,1:3)
    plot(xdata,ydata,'o',xdata,y0,'g',xdata,yfit,'b')
    text(0.02,0.90,txt_cen,'sc')
    text(0.02,0.85,txt_fwhm,'sc')
    text(0.02,0.80,txt_bg,'sc')
    subplot(4,1,4)
    plot(xdata,rs/pfit(1),'-x')
    waitforbuttonpress
    end
end
end
%%

% %% fit TOA by ceria peak
pk_use = [1:4];
xdata = dhkl(pk_use)';               % ideal d-spacing
list = 5:14;
ydata = p.cen(pos,list(pk_use))*(ChToEnergy(1))+ChToEnergy(2);   
yobserv = p.cen(pos,list(pk_use))*(ChToEnergy(1))+ChToEnergy(2);    % observed Energy (KeV)

fn = {'eddBragg2'}; npar = [1];
pfix = [nan]; % fix fitting parameters (currently none are fixed)
opt = optimoptions('lsqnonlin', 'Display','off','Jacobian','off','tolx',1e-8,'tolf',1e-8,'MaxIter',500);
p0  = [3.0];
pUB = [5.5];
pLB = [2.5];
y0   = sumfun1(p0,pfix,npar,fn,xdata);
            
       
[pfit,rn,rs,ex,out,lam,jac] = lsqnonlin('sumfun1',p0,pLB,pUB,opt,pfix,npar,fn,xdata,yobserv,ones(size(yobserv)));

ycalc  = sumfun1(pfit,pfix,npar,fn,dhkl(pk_use));   % ideal Energy (KeV)

figure(8);clf('reset');

plot(eddBragg2(pfit,dhkl(pk_use)),yobserv-ycalc','marker','o')
line(Ce_emission(1,1:3)/1000,[p.cen(pos,1:3)*ChToEnergy(1)+ChToEnergy(2)]-Ce_emission(1,1:3)/1000,'marker','^','color','g','linestyle','none','displayname','Ce-emission');
xlabel('Energy (KeV)');
ylabel('\DeltaE (KeV)');
title(sprintf('TOA = %2.15f, rn = %f',pfit(1),rn))
xlim([30 150])
ylim([-0.2 0.2])

text(0.2,0.95,sprintf('[%2.15f %2.15f]',ChToEnergy(1),ChToEnergy(2)),'sc');

TOA = pfit;

disp(sprintf('Use the following parameters:'));
disp(sprintf('Ch_to_Energy = [%2.15f %2.15f]',ChToEnergy(1),ChToEnergy(2)));
disp(sprintf('TOA = %2.15f (rnx1E3=%1.8f)',pfit(1),rn*1E3));


if if_plot(2)
    x = Ch*ChToEnergy(1)+ChToEnergy(2);
    figure(2); clf('reset')
    set(gcf,'Position',[200 400 800 600])
    line(x,sum(cegv.data{detno}(sum_range,:),1)/length(sum_range)/cegv.exp_time,'color',cc(3,:),'displayname','CeO2-GV-scan');   %
    line(Ehkl,ones(length(Ehkl),1),'marker','^','color','b','linestyle','none','displayname','CeO2\_hkl');
    line(Ce_emission(1,:)/1000,Ce_emission(2,:),'marker','^','color','g','linestyle','none','displayname','Ce-emission');
    
    set(gca,'yscale','log');
    grid('on')
    xlim([20 150])
    xlabel('Energy (KeV)');
    ylabel('Intensity (cps)');
    legend('toggle')
end
