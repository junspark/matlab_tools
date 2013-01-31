function [peakInfo]=fitModelDataP(polImg,modelPeakBounds,instr,material,cakeParms)
numPeaks=0;
for ii=1:length(material.latParm)
    numPeaks=numPeaks+length(material.tth{ii});
end

warning off MATLAB:polyfit:RepeatedPointsOrRescale
if size(polImg.radius,1)==1
    ttt=[];
    for ii=1:cakeParms.bins(1)
        ttt=[ttt; polImg.radius];
    end
    polImg.radius=ttt;
end

peakInfo.amp={};
peakInfo.mix={};
peakInfo.loc={};
peakInfo.G={};
StepSize=instr.StepSize;
WindowSize=instr.WindowSize;
Base=instr.Base;
Levels=instr.Levels;
FWHHFilter=instr.FWHHFilter;
HeightFilter=instr.HeightFilter;
PeakLocation=instr.PeakLocation;
for jj=1:cakeParms.bins(1)
    xdata=polImg.radius(jj,:);
    Akeep=[];
    lockeep=[];
    Mixkeep=[];
    Gkeep=[];
    
    ydata=msbackadj(xdata',polImg.intensity(jj,:)','StepSize',StepSize,'WindowSize',WindowSize,'ShowPlot','false','SmoothMethod','rloess');
    
    for ii=1:size(modelPeakBounds,1)
        [aa bb]=min(abs(polImg.radius(jj,:)-modelPeakBounds(ii,1)/instr.pixelSize));
        [cc dd]=min(abs(polImg.radius(jj,:)-modelPeakBounds(ii,2)/instr.pixelSize));
        [A Aind]=max(ydata(bb:dd));
        rng=bb:dd;
        
        opts = optimset(...
            'Display', 'off', ...
            'DerivativeCheck', 'off', ...
            'TolFun', .1e-14, ...
            'TolX', .8e-15, ...
            'MaxIter', 1e3,...
            'MaxFunEvals',2e5, ...
            'Display','off' ...
            );
        
        x0=[A .4 6 xdata(rng(Aind)) 0 0];
        lb=[0 0 0 xdata(rng(Aind))-20 -inf -inf];
        ub=[A*10 1 10 xdata(rng(Aind))+20 inf inf];
        %         lb=[0 0 0 xdata(rng(1)) -inf -inf];         % changed by Atish
        %         ub=[A*10 1 10 xdata(rng(end)) inf inf];     % changed by Atish
        xOut=lsqcurvefit(@peakRefinementP,x0,xdata(rng),ydata(rng)',lb,ub,opts);
        
        Akeep=[Akeep xOut(1)];      % peak intensity or amplitude
        Mixkeep=[Mixkeep xOut(2)];  % mixing parameter
        Gkeep=[Gkeep xOut(3)];      % FWHM (in pixels)
        lockeep=[lockeep xOut(4)];  % location for the refined peak (in pixels)
        
        figure(100)
        clf
        plot(xdata(rng), ydata(rng)', 'k.')
        hold on
        plot(xdata(rng), peakRefinementP(x0, xdata(rng)), 'r-')
        plot(xdata(rng), peakRefinementP(xOut, xdata(rng)), 'g-')
        hold off
        xlabel('radial distance (pixels)')
        ylabel('intensity (cts)')
        pause(0.1)
    end
    peakInfo.amp{jj}=Akeep;
    peakInfo.mix{jj}=Mixkeep;
    peakInfo.G{jj}=Gkeep;
    peakInfo.loc{jj}=lockeep;
end


peakInfo.amp=cell2mat(peakInfo.amp');
peakInfo.mix=cell2mat(peakInfo.mix');
peakInfo.G=cell2mat(peakInfo.G');
peakInfo.loc=cell2mat(peakInfo.loc');

