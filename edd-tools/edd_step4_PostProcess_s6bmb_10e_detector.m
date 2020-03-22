clear all
close all
clc

%%% PEAK TO FIT
d_pk    = [3.45 2.80];
d_UB    = [3.65 3.04];
d_LB    = [3.21 2.54];

edf_cal = Read_s6bmb_edf('./ALU_021820_0001.edf');
for i = 1:1:10
    b(i)    = edf_cal(i).polynomial(1);
    m(i)    = edf_cal(i).polynomial(2);
    tth(i)  = edf_cal(i).tth(1);
end
for mmm = 1:1:length(d_pk)
    Epk(mmm,:)  = Angstrom2keV(2*d_pk(mmm)*sind(tth/2));
end

proot   = './data';
% froot   = 'b949_81_84_longscan000';
froot   = 'm592_62_longscan000';

pfroot  = fullfile(proot, froot, sprintf('%s.xy', froot));
edd     = readedd5_6bm(pfroot);

numscans    = length(edd);

nx  = 16;
ny  = 10; 
nz  = 17;
for iii = 1:1:length(d_pk)
    samXE   = [];
    samYE   = [];
    samZE   = [];
    
    Aout    = [];
    BkgInt  = [];
    Eout    = [];
    Gout    = [];
    PkInt   = [];
    ef      = [];
    nout    = [];
    rsn     = [];
    
    %%% 160 for b949_81_84_longscan000
    %%% 100 for m592_62_longscan000
    for jjj = 1:1:100 %numscans 
        fname_fits  = sprintf('%s_scan%d_fits.mat', froot, jjj);
        
        fitdata = load(fname_fits);
        
        samXE   = [samXE; fitdata.samXE];
        samYE   = [samYE; fitdata.samYE];
        samZE   = [samZE; fitdata.samZE];
        
        Aout    = [Aout; cell2mat(fitdata.Aout(:,iii))'];
        BkgInt  = [BkgInt; cell2mat(fitdata.BkgInt(:,iii))'];
        Eout    = [Eout; cell2mat(fitdata.Eout(:,iii))'];
        Gout    = [Gout; cell2mat(fitdata.Gout(:,iii))'];
        PkInt   = [PkInt; cell2mat(fitdata.PkInt(:,iii))'];
        ef      = [ef; cell2mat(fitdata.ef(:,iii))'];
        nout    = [nout; cell2mat(fitdata.nout(:,iii))'];
        rsn     = [rsn; cell2mat(fitdata.rsn(:,iii))'];
    end
    
    [numpts, numdet]    = size(Aout);
    
    plotItem.data       = PkInt-BkgInt;
    plotItem.name       = 'intensity';
    % plotItem.datarange  = [min(plotItem.data(:)) max(plotItem.data(:))];
    % plotItem.datarange  = [0 100; 0 400]; %%% for b949_81_84_longscan000
    plotItem.datarange  = [0 150; 0 1200]; %%% for m592_62_longscan000
    
%     plotItem.data       = rsn;
%     plotItem.name       = 'ef';
%     plotItem.datarange  = [min(plotItem.data(:)) 5000];
    
%     plotItem.data       = ef;
%     plotItem.name       = 'ef';
%     plotItem.datarange  = [min(plotItem.data(:)) 5];
    
%     plotItem.data       = repmat(Epk(iii,:), numpts, 1)./Eout - 1;
%     plotItem.name      = 'strain';
%     plotItem.datarange  = [-1e-2 1e-2];
    
    % idx_plot    = ef > 0;
    idx_plot    = ef < 4;       %%% THIS PLOTS ALL
    for jjj = 1:1:numdet
        idx_jjj = idx_plot(:,jjj);
        
        fname_gcf1  = sprintf('%s_%s_pk%d_detno%d_dots.png', froot, plotItem.name, iii, jjj);
        fname_gcf2  = sprintf('%s_%s_pk%d_detno%d_contour.png', froot, plotItem.name, iii, jjj);
        fname_csv   = sprintf('%s_%s_pk%d_detno%d.csv', froot, plotItem.name, iii, jjj);
        
        %%%%
        figure(1)
        scatter3(samXE(idx_jjj), samYE(idx_jjj), samZE(idx_jjj), ...
            [], plotItem.data(idx_jjj,jjj), 'filled')
        colorbar vert
        colormap jet
        axis equal tight
        caxis(plotItem.datarange(iii, :))
        title(sprintf('%s | det_num=%d | using pk at %2.2f Ang', plotItem.name, jjj, d_pk(iii)), 'Interpreter', 'none')
        xlabel('samX (mm)')
        ylabel('samY (mm)')
        zlabel('samZ (mm)')
        view([35 -36])
        saveas(gcf, fname_gcf1, 'png')
        
        %%%%
        F   = scatteredInterpolant(samXE(idx_jjj), samYE(idx_jjj), samZE(idx_jjj), plotItem.data(idx_jjj,jjj), ...
            'natural');
        
        xxx = min(samXE(idx_jjj)):0.5:max(samXE(idx_jjj));
        yyy = min(samYE(idx_jjj)):0.5:max(samYE(idx_jjj));
        zzz = min(samZE(idx_jjj)):0.5:max(samZE(idx_jjj));
        
        [xq, yq, zq]    = meshgrid(xxx, yyy, zzz);
        cq              = F(xq,yq,zq);
        
        xslice = mean(xxx); 
        yslice = mean(yyy); 
        zslice = mean(zzz);
        
        vm  = [35 -36; ...
            0 90; ...
            0 0; ...
            90 0; ...
            ];
        
        figure(2)
        set(gcf, 'Position', [0 10 1000 1000])
        for kkk = 1:1:4
            subplot(2, 2, kkk)
            slice(xq,yq,zq,cq,xslice,yslice,zslice)
            colorbar vert
            colormap jet
            axis equal tight
            caxis(plotItem.datarange(iii, :))
            title(sprintf('%s | det_num=%d | using pk at %2.2f Ang', plotItem.name, jjj, d_pk(iii)), 'Interpreter', 'none')
            xlabel('samX (mm)')
            ylabel('samY (mm)')
            zlabel('samZ (mm)')
            view(vm(kkk, :))
        end
        saveas(gcf, fname_gcf2, 'png')
        
        close all
        %%%%
        Aoutjjj     = Aout(:,jjj);
        BkgIntjjj   = BkgInt(:,jjj);
        Eoutjjj     = Eout(:,jjj);
        Goutjjj     = Gout(:,jjj);
        PkIntjjj    = PkInt(:,jjj);
        PkIntnetjjj = plotItem.data(:,jjj);
        efjjj       = ef(:,jjj);
        noutjjj     = nout(:,jjj);
        rsnjjj      = rsn(:,jjj);
        
        tableout    = table(samXE, samYE, samZE, ...
            Aoutjjj, BkgIntjjj, Eoutjjj, Goutjjj, PkIntjjj, PkIntnetjjj, ...
            efjjj, noutjjj, rsnjjj);
        
        writetable(tableout, fname_csv)
    end
%     pause
end
return
