% eddplot to visdualize 6BM data.
%
% Copyright 2017 Andrew Chuang (chuang.cp@gmail.com)
% $Revision: 1.1 $  $Date: 2017/09/20 $

function eddplot2(da,opt)

fsa = 13;
fst = 18;
cc = lines(20);

if nargin == 1
    opt = '';
end

%%%%%%  define default option %%%%%

if isfield(opt,'datalim')
    datalim = opt.datalim;
else
    % range for (Ch/E/d),(posno),(Intensity)
    datalim = {'auto','auto','auto'}; 
end

if isfield(opt,'normalize')
    normalize = opt.normalize;
else
    normalize = 0;
end

if isfield(opt,'detno')
    detno = opt.detno;
else
    detno = 1;
end

if isfield(opt,'phase')
    phase = opt.phase;
else
    phase = 1;
end

if isfield(opt,'pk')
    pk = opt.pk;
else
    pk = 5;
end

if isfield(opt,'do_export')
    do_export = opt.do_export;
else
    do_export = 0;
end

if isfield(opt,'x_unit')
    x_unit = lower(opt.x_unit);
else
    x_unit = 'ch';
end

if isfield(opt,'x_range')
    x_range = opt.x_range;
else
    x_range = ':';
end

% if isfield(opt,'fig_handle')&&ishandle(opt.fig_handle)
%     fig_handle = opt.fig_handle;
% else
%     fig_handle = 0;
% end

if isfield(opt,'avg_1D')
    avg_1D = opt.avg_1D;
else
    avg_1D = 0;
end

if isfield(opt,'title')
    title_text = opt.title;
else
    title_text = '';
end

if isfield(opt,'dzero')
    dzero = opt.dzero;
else
    dzero = zeros(1,length(pk))+1;
end
if isfield(opt,'posno')
    posno = opt.posno;
else
    posno = 1;
end

if isfield(opt,'xoffset')
    xoffset = opt.xoffset;
else
    xoffset = 0;
end

if isfield(opt,'type')
    type = opt.type;
else
    type = '2draw';
end

if isfield(opt,'scno')
    scno = opt.scno(1);   % only plot one scan at a time. (May.17)
else
    scno = 1;
end

if isfield(da,'Material')
    label_peak = 1;
else
    label_peak = 0;
end

if ~strcmp(type,'raw')&isfield(da(1),'Material')
    pkname = da(1).Material(phase).hkls;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type,'2draw')&&size(da(scno).data{detno},1)==1
    type = 'raw';
end

if isfield(opt,'Inst')
    da(1).Inst = opt.Inst;
elseif isfield(da(scno),'Inst') && length(da(scno).Inst)==1
    da(scno).Inst(2).Ch2E = da(1).Inst(1).Ch2E;
    da(scno).Inst(2).TOA  = da(1).Inst(1).TOA; 
elseif ~isfield(da(scno),'Inst')
    fprintf('No Instrument parameters provided!!\n');
    fprintf('Use default TOA  = 3 (deg)\n');
    fprintf('Use default Ch2E = [0.03 0];\n');
    da(scno).Inst(1).Ch2E = [0.03 0];
    da(scno).Inst(1).TOA  = 3;
    da(scno).Inst(2).Ch2E = [0.03 0];
    da(scno).Inst(2).TOA  = 3;
end

% calculate E_grid, d_grid
hc = 12.398419057638671;
for i = 1:2
    da(scno).Inst(i).E_grid = [1:8192]*da(1).Inst(i).Ch2E(1)+da(1).Inst(i).Ch2E(2);
    da(scno).Inst(i).d_grid = hc./da(scno).Inst(i).E_grid*0.5/sind(da(1).Inst(i).TOA/2);
end

%%%%%% switch to det-2 data
switch detno
    case 1
        if strcmp(type,'2draw')
            hfig = findall(0,'Tag','edd_fig_det1_map');
            if ishandle(hfig)
                fig = hfig(1); 
                clf(fig,'reset');
                set(fig,'Tag','edd_fig_det1_map');
            else
                fig = figure(161);
                set(fig,'Position',[50 100 800 600],'Tag','edd_fig_det1_map');
            end
            figure(161);
        else
            hfig = findall(0,'Tag','edd_fig_det1_line');
            if ishandle(hfig)
                fig = hfig(1); 
                clf(fig,'reset');
                set(fig,'Tag','edd_fig_det1_line');
            else
                fig = figure(162);
                set(fig,'Position',[50 700 800 600],'Tag','edd_fig_det1_line');
            end
            figure(162);
        end
    case 2
        if strcmp(type,'2draw')
            hfig = findall(0,'Tag','edd_fig_det2_map');
            if ishandle(hfig)
                fig = hfig(1); 
                clf(fig,'reset');
                set(fig,'Tag','edd_fig_det2_map');
            else
                fig = figure(261);
                set(fig,'Position',[100 100 800 600],'Tag','edd_fig_det2_map');
            end
            figure(261);
        else
            hfig = findall(0,'Tag','edd_fig_det2_line');
            if ishandle(hfig)
                fig = hfig(1); 
                clf(fig,'reset');
                set(fig,'Tag','edd_fig_det2_line');
            else
                fig = figure(262);
                set(fig,'Position',[100 700 800 600],'Tag','edd_fig_det2_line');
            end
            figure(262);
        end        
%         for i = 1:length(da)
%             da(i).data = da(i).data2;
%         end

        %da(scno).data = cat(3,da(scno).data,da(scno).data2);
    otherwise
        fprintf('3 or more detector config is not supported yet!!\n');
        return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start ploting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc = lines(100);
ax = axes('parent',hfig,'fontsize',fsa,'box','on');
grid on;

switch lower(type)
    case 'cen'
        xdata = da(scno).motorpos(x_range)+xoffset;
        display(sprintf('\nAveraged Peak Position (Mean%sStd) for %s',char(177),title_text));
        for i = 1:length(pk)
            ydata = da(scno).fit(phase,pk(i)).cen(x_range);
            flag  = find(ydata~=0);
            ydata = ydata(flag);
            if normalize == 1
                line(xdata(flag),((mean(ydata)./ydata)-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{-6})');
                %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
                display(sprintf('[%s] = %8.2f %s %8.2f, ',num2str(pkname(pk(i),:)),mean(((ydata./mean(ydata))-1)*1E6),char(177),std(((ydata./mean(ydata))-1)*1E6)));
            elseif normalize == 0
                line(xdata(flag),ydata,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Peak Center (keV)');
                display(sprintf('[%s] = %8.8f %s %6.4f, ',num2str(pkname(pk(i),:)),mean(ydata),char(177),std(ydata)));
            else
                line(xdata(flag),((dzero(pk(i))./ydata)-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{-6})');
                display(sprintf('[%s] = %8.2f %s %8.2f, ',num2str(pkname(pk(i),:)),mean((ydata./dzero(pk(i))-1).*1E6),char(177),std((ydata./dzero(pk(i))-1).*1E6)));
                %std((ydata./dzero(pk(i))-1).*1E6)
            end
            %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
            
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1});
        ylim(datalim{3});
        title(title_text,'fontsize',fst);
        leg = legend('toggle');

    case 'dspac'
        xdata = da(scno).motorpos(x_range)+xoffset;
        display(sprintf('\nAveraged Peak Position (Mean%sStd) for %s',char(177),title_text));
        for i = 1:length(pk)
            %da(scno).fit(phase,pk(i),detno);
            ydata = da(scno).fit(phase,pk(i),detno).dspac(x_range);
            flag  = find(ydata~=0);
            ydata = ydata(flag);
            if normalize == 1
                line(xdata(flag),((mean(ydata)./ydata)-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{-6})');
                %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
                display(sprintf('[%s] = %8.2f %s %8.2f, ',num2str(pkname(pk(i),:)),mean(((ydata./mean(ydata))-1)*1E6),char(177),std(((ydata./mean(ydata))-1)*1E6)));
            elseif normalize == 0
                line(xdata(flag),ydata,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('d-spacing (Angstron)');
                display(sprintf('[%s] = %8.8f %s %6.4f, ',num2str(pkname(pk(i),:)),mean(ydata),char(177),std(ydata)));
            else
                line(xdata(flag),((ydata./dzero(pk(i)))-1)*1E6,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                %datalim{2} = [-5E3 5E3];
                ylabel('Strain (x 10^{6})');
                display(sprintf('[%s] = %8.2f %s %8.2f, ',num2str(pkname(pk(i),:)),mean((ydata./dzero(pk(i))-1).*1E6),char(177),std((ydata./dzero(pk(i))-1).*1E6)));
                %std((ydata./dzero(pk(i))-1).*1E6)
            end
            %display(sprintf('[%s] = %04.8f, ',num2str(pkname(pk(i),:)),mean(ydata)));
            
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1});
        ylim(datalim{3});
        title(title_text,'fontsize',fst);
        leg = legend('toggle');        
        
    case 'int'
        xdata = da(scno).motorpos+xoffset;
        for i = 1:length(pk)
            ydata = da(scno).fit(phase,pk(i)).int;
            flag  = find(ydata~=0);
            ydata = ydata(flag);
            if normalize == 0
                line(xdata(flag),ydata.*da(1).exp_time,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Intensity (counts)');
            else
                line(xdata(flag),ydata,'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
                ylabel('Normalized Intensity (cps)');
            end
            text(0.4,0.95,sprintf('Total count time: %d sec',da(1).exp_time),'sc','fontsize',fsa);
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1})
        ylim(datalim{3})
        title(title_text,'fontsize',fst)
        leg = legend('toggle');        
        
    case 'fwhm'
        xdata = da(scno).motorpos+xoffset;
        for i = 1:length(pk)
            ydata = da(scno).fit(phase,pk(i)).fwhm;
            flag  = find(ydata~=0);
            line(xdata(flag),ydata(flag),'marker','o','color',cc(i,:),'displayname',sprintf('%s',num2str(pkname(pk(i),:))));
            ylabel('FWHM (keV)');
        end
        set(gca,'fontsize',fsa);
        xlabel('Position (mm)');
        xlim(datalim{1})
        ylim(datalim{3})
        title(title_text,'fontsize',fst)
        leg = legend('toggle');        
        
    case 'raw'
        yrange = datalim{3};
        switch x_unit
            case 'd'
                %xxdata = da(point).Inst.d_grid;
                xxdata = da(scno).Inst(detno).d_grid;
                xrange = [0 10];
                xlab = sprintf('d (%s)',char(197));
            case 'e'
                %xxdata = da(point).Inst.E_grid;
                xxdata = da(scno).Inst(detno).E_grid;
                xrange = [0 200];
                xlab = sprintf('E (keV)');
            case 'ch'
                xxdata = 1:8192;
                xrange = [1 8192];
                xlab = sprintf('Channel');
        end
        if avg_1D
            if normalize == 0
                yydata = sum(da(scno).data{detno}(posno,:),1)./length(posno);
                tx_yax = 'Intensity (counts)';
            else
                yydata = sum(da(scno).data{detno}(posno,:),1)./length(posno)./da(scno).exp_time;
                tx_yax = 'Normalized Intensity (cps)';
            end
            line(xxdata, yydata,'color',cc(1,:),'marker','.','displayname','raw curve')
        else
            if normalize == 0
                yydata = da(scno).data{detno}(posno,:);
                tx_yax = 'Intensity (counts)';
            else
                yydata = da(scno).data{detno}(posno,:)./da(scno).exp_time;
                %yydata = sum(da(scno).data(posno,:),1)./length(posno);
                tx_yax = 'Normalized Intensity (cps)';
            end
            for k = 1:length(posno)
                line(xxdata, yydata(k,:),'color',cc(k,:),'marker','.','displayname',sprintf('pos-%d',posno(k)))
            end
        end

%         if strcmp(x_unit,'d')
%             xxx = 0.5:0.001:10;
%             yyy = interp1(xxdata,yydata,xxx,'spline');
%             line(xxx,yyy,'color','r','marker','o','markersize',3);
%         end
        
%         if label_peak
%             ccode = {'b' '[0 0.8 0]' '[0.73 0.7 0.7]'};
%             %pk_list = [da(point).pkid_fit{:}];
%             pk_list = [da(1).pkid_fit{:}];
%             for m = 1:length(da(1).Material)
%                 phase_name_list{m} = da(1).Material(m).name;
%                 switch x_unit
%                     case 'd'
%                         xxsample = da(1).Material(m).d_hkl;
%                     case 'E'
%                         xxsample = da(1).Material(m).E_hkl;
%                     case 'Ch'
%                         xxsample = da(1).Material(m).Ch_hkl;
%                 end
%                 line(xxsample, ones(length(xxsample), 1),'marker','^','color',ccode{m},'markersize',9,'linestyle','none','Displayname',sprintf('%s',da(1).Material(m).name))
%                 lnfit = line(xxsample(pk_list([da(1).phase_fit{:}]==m)), ones(length(pk_list([da(1).phase_fit{:}]==m)), 1),'marker','o','color','r','markersize',5,'linestyle','none','displayname','peak fitted');
%                 if m ~= 1; set(get(get(lnfit,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');end
%             end
%         end
        if ~strcmp(datalim{1},'auto'); xrange = datalim{1};end
        grid on
        set(gca,'yscale','log','fontsize',fsa,'box','on');
        xlabel(xlab)
        ylabel(tx_yax)
        xlim(xrange);
        ylim(yrange);
        legend('toggle');
        if length(posno)>1
            title_text = sprintf('Det-%d, Positon = %d - %d',detno,posno(1),posno(end));
        else
            title_text = sprintf('Det-%d, Positon = %d',detno,posno(1));
        end
        title(title_text,'fontsize',fst)
        text(0.1,0.05,sprintf('Total count time: %d sec',da(1).exp_time*length(posno)),'sc','fontsize',fsa);
        
    case '2draw'
        switch x_unit
            case 'ch'
                xrange = datalim{1};
                yrange = datalim{2};  %yrange = 'auto';
                int_range = datalim{3};
                if ~strcmp(yrange,'auto'); yrange=yrange+[-0.5 0.5];end
                xlab = sprintf('Channel');
                xgrid =1:8192;
                ygrid =1:size(da(scno).data{detno}(:,:),1);
            case 'e'
                xrange = datalim{1};
                yrange = datalim{2};  %yrange = 'auto';
                int_range = datalim{3};
                if ~strcmp(yrange,'auto'); yrange=yrange+[-0.5 0.5];end
                xlab = sprintf('Energy (keV)');
                xgrid= da(scno).Inst(detno).E_grid;
                %xgrid = [1:8192]*Inst(detno).Ch2E(1)+Inst(detno).Ch2E(2);
                ygrid = 1:size(da(scno).data{detno}(:,:),1);
            case 'd'
                xrange = datalim{1};
                yrange = datalim{2};
                int_range = datalim{3};
                if ~strcmp(yrange,'auto'); yrange=yrange+[-0.5 0.5];end
                xlab = sprintf('d (Angstron)');
                if ~isfield(da(scno),'d_data')
                    for k = 1:size(da(scno).data{detno}(:,:),1)
                        da(scno).Inst(detno).d_grid2 = 0.4:0.001:10;
                        xx = da(scno).Inst(detno).d_grid;
                        yy = da(scno).data{detno}(k,:);
                        qy = interp1(xx,yy,da(scno).Inst(detno).d_grid2,'spline');
                        qy(qy<0.1|qy>max(yy)*2)=0;
                        da(scno).d_data(k,:,detno) = qy;
                    end
                end
                xgrid= da(scno).Inst(detno).d_grid2;
                ygrid = 1:size(da(scno).d_data(:,:,detno),1);
                %assignin('base','da',da);
            otherwise
                fprintf('The desire x_unit: %s is not supported yet!!\n',x_unit);
                return;
        end
        ylab = sprintf('%s (#) steps',da(scno).motorname);
        linesap = find(da(scno).log==char(10));
        wordsap = find(da(scno).log(1:linesap(1)-1)==' ');
        title_text=sprintf('%s',da(scno).log(wordsap(3)+1:linesap(1)-1));
        
        if strcmp(x_unit,'d')
            imagesc(xgrid,ygrid,log10(da(scno).d_data(:,:,detno)));
            %set(gca,'clim',[]);
        else
            imagesc(xgrid,ygrid,log10(da(scno).data{detno}(:,:)));
        end
        
        if ~strcmp(xrange,'auto'); xlim(xrange);end
        if ~strcmp(yrange,'auto'); ylim(yrange);end
        if ~strcmp(int_range,'auto'); set(gca,'clim',[log10(int_range(1)) log10(int_range(2))]);end
        xlabel(xlab,'fontsize',13);
        ylabel(ylab,'fontsize',13);
        title(title_text,'fontsize',17);
    otherwise
        fprintf('Selected plot type is not supported');
        
end