if rebinImg
    fileStps    = length(darkNum);
    for i = 1:1:fileStps
        file.name   = darkList{i};
        [a b c]     = fileparts(file.name);
        instr.fileName  = [file.path b];
        if i == 1
            darkData    = imread([instr.fileName '.tiff']);
        else
            darkData    = darkData + imread([instr.fileName '.tiff']);
        end
    end
    darkData    = darkData./fileStps;
    
    % FOR CHESS GE
    % darkData    = imrotate(darkData, -90);
    darkData    = fliplr(darkData);
end

fileStps    = length(fileNum);
for i = 1:1:fileStps
    file.name   = fileList{i};
    [a b c]     = fileparts(file.name);
    instr.fileName  = [file.path b];
    
    if makeTIFs
        exeMARCVT   = ['!./Mar345/marcvt -tiff16 -colors 65536 -min 0 -max 65536 ' file.path file.name];
        eval(exeMARCVT)
        imageData=imread([instr.fileName '.tiff16']);
        delImg=['!rm ' instr.fileName '.tiff16'];
        eval(delImg)
    end
    
    if rebinImg
        [a b c]         = fileparts(file.name);
        instr.fileName  = [file.path b];
        
        imageData	= imread([instr.fileName '.tiff']);
        
        % FOR CHESS GE
        % imageData   = imrotate(imageData, 90);
        imageData   = fliplr(imageData);
        
        imageData   = imageData - darkData;
        
%         figure(101)
%         clf
%         imagesc(log(double(imageData)))
%         axis square tight
    end
    
    if useCeO2instr
        load([file.path CeO2List{i}])
        instr=par.instr;
        instr.fileName  = [file.path b];
        instr.distance  = instr.distance-1;
    else
        instr.omega = stndOmega(i);
        instr.chi   = stndChi(i);
    end
    
    if rebinImg
        if i == 1
            mesh    = buildMeshXRD(imageData);
        end
        polImg  = reBin(mesh,instr,cakeParms,imageData);
        
        save([file.path b '.polImg.mat'],'polImg');
    else
        instr.fileName=[file.path b];
        disp([file.path b '.polImg.mat']);
        load([file.path b '.polImg.mat']);
    end
    [material,hkls] = findModelPeaks(instr,material,polImg);
    peakInfo        = fitModelDataP(polImg,material.rList,instr,material,cakeParms);
    pause(2)
    
    %0 no cor, 1 mar, 2 GE, 3 GE for the mar
    if instr.DistType == 0
        DistPars=[];
    elseif instr.DistType == 1
        DistPars=instr.a0;
    elseif instr.DistType == 2
        DistPars=[instr.a1 instr.a2 instr.n1 instr.n2 instr.t1 instr.e1 instr.detectorSize];
    elseif instr.DistType == 3
        DistPars=[instr.a1 instr.a2 instr.n1 instr.n2 instr.t1 instr.e1 instr.detectorSize];
    end
    
    opts = optimset(...
        'DerivativeCheck', 'off', ...
        'TolFun', 1e-18, ...
        'TolX', 1e-22, ...
        'MaxIter', 10000,...
        'MaxFunEvals',3e5,...
        'Display','final',...
        'TypicalX',[10 -10 1 1 1 DistPars]...
        );
    
    assignin('base','instr',instr);
    assignin('base','peakInfo',peakInfo);
    assignin('base','material',material);
    
    x0      = [instr.centers(1) instr.centers(2) instr.gammaY instr.gammaX 1 DistPars];
    funcOut = peakOnlyGeomModel(x0);
    [x,Resnorm,FVAL,EXITFLAG,OUTPUT]    = lsqnonlin(@peakOnlyGeomModel,x0,[],[],opts);%
    
    instr.centers   = [x(1)    x(2)];
    instr.gammaY    = x(3);
    instr.gammaX    = x(4);
    instr.distance  = x(5)*instr.distance;
    DistPars        = x(6:end);
    
    % 1 mar, 2 GE, 3 GE for the mar
    if instr.DistType == 1
        instr.a0=DistPars(1);
        DistPars=instr.a0;
    elseif instr.DistType == 2
        instr.a1=DistPars(1);
        instr.a2=DistPars(2);
        instr.n1=DistPars(3);
        instr.n2=DistPars(4);
        DistPars=[instr.a1 instr.a2 instr.n1 instr.n2 instr.t1 instr.e1 instr.detectorSize];
    elseif instr.DistType == 3
        instr.a1=DistPars(1);
        instr.a2=DistPars(2);
        instr.n1=DistPars(3);
        instr.n2=DistPars(4);
        instr.t1=DistPars(5);
        instr.e1=DistPars(6);
        DistPars=[instr.a1 instr.a2 instr.n1 instr.n2 instr.t1 instr.e1 instr.detectorSize];
    end
    
    mappedData = XFormRadialSpectraPeaks(instr.centers./1000,instr.distance, ...
        instr.gammaY/1000, instr.gammaX/1000, ...
        peakInfo.loc'.*instr.pixelSize, ...
        cakeParms.azim, instr.tthRng,...
        DistPars,instr.DistType);
    
    %%%%make plots to show the quality of the results
    qual=mappedData';
    rmapped=instr.distance.*tand(mappedData);
    rA=instr.distance.*tand(repmat(material.tth{1}(material.modelPeaks),cakeParms.bins(1),1))';
    figure(2);
    clf;
    tmp=1.96*std(sind(repmat(material.tth{1}(material.modelPeaks),cakeParms.bins(1),1)./2)./sind(qual./2)-1);
    [AX,H1,H2] = plotyy(material.tth{1}(material.modelPeaks),tmp,material.tth{1}(material.modelPeaks),mean(peakInfo.amp),'plot');
    set(get(AX(1),'Ylabel'),'String','U_i^m');
    set(get(AX(2),'Ylabel'),'String','Amp.');
    xlabel('2\theta');
    set(H1,'Marker','.');
    set(H2,'Marker','s');
    figure(3);clf;plot(material.tth{1}(material.modelPeaks),sind(repmat(material.tth{1}(material.modelPeaks),cakeParms.bins(1),1)./2)./sind(qual./2)-1,'.');
    ylabel('Pseudo Lattice Strain','FontSize',14,'FontWeight','Bold');
    xlabel('2\theta (deg.)','FontSize',14,'FontWeight','Bold');
    A=gca;
    get(A);
    set(A,'FontSize',14,'FontWeight','Bold');
    legend(num2str(hkls(1:length(material.tth{1}),:)));
    figure(4);clf;plot(1:cakeParms.bins(1),(rmapped-rA)./rA,'.');
    ylabel('\rho-\rho_e_x_p_e_c_t_e_d','FontSize',14,'FontWeight','Bold');
    xlabel('\eta bin number','FontSize',14,'FontWeight','Bold');
    A=gca;
    get(A);
    set(A,'FontSize',14,'FontWeight','Bold');
    legend(num2str(hkls(1:length(material.tth{1}),:)));
    
    aa      = [];
    azim    = cakeParms.azim;
    for  j = 1:length(azim)
        rdata   = tand(mappedData(:,j))*instr.distance;
        
        if instr.DistType==0
            rhoHatdPrime=rdata;
        elseif instr.DistType==1
            rhoHatdPrime  = rdata + DistPars(1);
        elseif instr.DistType==2
            f = DistPars(1)*((rdata/DistPars(7)).^(DistPars(3)))*cosd(DistPars(5)*azim(i)+DistPars(6))...
                + DistPars(2)*((rdata/DistPars(7)).^DistPars(4)) + 1;
            rhoHatdPrime  = f.*rdata;
        elseif instr.DistType==3
            f = DistPars(1)*((rdata/DistPars(7)).^(DistPars(3)))*cosd(DistPars(5)*azim(i)+DistPars(6))...
                + DistPars(2)*((rdata/DistPars(7)).^DistPars(4)) + 1;
            rhoHatdPrime  = f.*rdata;
        end
        aa  = [aa rhoHatdPrime];
    end
    
    figure(5);clf;plot(1:cakeParms.bins(1),aa-repmat(tand(material.tth{1}(material.modelPeaks)')*instr.distance,1,length(azim)),'.')
    legend(num2str(hkls(1:length(material.tth{1}),:)))
    ylabel('distortion cor. (mm)','FontSize',14,'FontWeight','Bold');
    xlabel('\eta bin number','FontSize',14,'FontWeight','Bold');
    A=gca;
    get(A);
    set(A,'FontSize',14,'FontWeight','Bold');
    legend(num2str(hkls(1:length(material.tth{1}),:)),'Location','NorthEastOutside');
    
    DATA=cell(1,cakeParms.bins(1));
    for ii=1:cakeParms.bins(1)
        DATA{ii}= [instr.pixelSize*polImg.radius(1,:)' polImg.intensity(ii,:)' ];
    end
    
    mappedData = XFormRadialSpectra(instr.centers/1000, instr.distance, ...
        instr.gammaY/1000, instr.gammaX/1000, ...
        DATA,cakeParms.azim, instr.tthRng,...
        DistPars,instr.DistType);
    
    %%%truncate the data to have a regular grid
    minMapped=1e10;
    maxMapped=0;
    for ii=1:cakeParms.bins(1)
        minMap=size(mappedData{ii},1);
        if minMap<minMapped
            minMapped=minMap;
        end
        maxMap=size(mappedData{ii},1);
        if maxMap>maxMapped
            maxMapped=maxMap;
        end
    end
    if maxMapped~=minMapped;
        'truncated edge of detector'
    end
    
    ttt=[];
    rrr=[];
    for ii=1:cakeParms.bins(1)
        rrr=[rrr; mappedData{ii}(1:minMapped,1)'];
        ttt=[ttt; mappedData{ii}(1:minMapped,2)'];
    end
    r=instr.distance.*tand(rrr);
    I=ttt;
    
    polImgCor.radius=r;
    polImgCor.intensity=I;
    polImgCor.tth=atand(polImgCor.radius./instr.distance);
    polImgCor.azimuth=polImg.azimuth;
    
    if generateESG
        buildESG(instr, cakeParms, polImgCor)
    end
    
    par.instr       = instr;
    par.material    = material;
    par.cakeParms   = cakeParms;
    par.peakInfo    = peakInfo;
    par.polImgCor   = polImgCor;
    par.qual        = qual;
    
    save([instr.fileName '.par.mat'],'par');
end