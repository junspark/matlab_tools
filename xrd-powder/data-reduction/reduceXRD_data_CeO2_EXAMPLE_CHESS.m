clear all
close all
clc

%%% DATA REDUCTION PARAMETERS
rebinImg        = 1;
refitPeaks      = 1;
generateESG     = 1;
GE              = 1;
useCeO2instr    = 0;

%%% INSTRUMENT PARAMETERS
instr.energy        = 59.390;                       % keV
instr.wavelength    = keV2Angstrom(instr.energy);   % wavelength (Angstrom)
instr.distance      = 941.4603;                     % mm
instr.centers       = [0 0];            % delta center x & y; USE [0 0] (NOT IMPORTANT UNLESS FOR STRAIN)
instr.gammaX        = 0;                % rad  %%% IF FIT2D ROTATIONS < 0 THEN USE 0 FOR BOTH (NOT IMPORTANT UNLESS FOR STRAIN)
instr.gammaY        = 0;                % rad
instr.radOffset     = 0;
instr.rock          = 5;
instr.tthRng        = [3 14];
instr.detectorSize  = 409.6;    %mm
instr.pixelSize     = 0.200;    %mm
instr.StepSize      = .5*100;
instr.WindowSize    = 100;
instr.Base          = 2;
instr.Levels        = 12;
instr.FWHHFilter    = .04;
instr.HeightFilter  = 900;
instr.PeakLocation  = 0.5;

instr.DistType      = 2 ;       % 0 no cor, 1 mar, 2 GE, 3 mar GE

%%% DETECTOR DISTORTION PARAMETERS - GE
instr.a0   = 0.0204;
instr.a1   = -0.0000504131769;
instr.a2   = -0.0010103284296;
instr.n1   = 1.3857;
instr.n2   = 2.1433;
instr.t1   = 2.0747;
instr.e1   = 17.9846;

%%% FILE DESIGNATIONS
file.path   = './example/CHESS/';  %% FILE PATH : NOTE "/" AT THE END
matBaseName = 'Ti_1MA_S2';

%%%
CeO2BaseName    = 'Ti_1MA_S2';
fileNum         = {...
    '_01113';'_01114'};
stndOmega   = [0 30];   % ROTATION ABOUT YL
stndChi     = [0 0];    % ROTATION ABOUT XL-PRIME
darkNum     = {'_01122';'_01123'};

%%% PAR.MAT LIST
CeO2List    = {};
for ii=1:length(stndOmega)
    CeO2List{ii}    = [CeO2BaseName fileNum{ii} '.par.mat'];
end

%%% FILE LIST
ct          = 1;
fileList    = {};
for ii = 1:length(fileNum)
    fileList{ct}    = [matBaseName fileNum{ii} '.tiff'];
    ct  = ct+1;
end

ct          = 1;
darkList    = {};
for ii=1:length(darkNum)
    darkList{ii}    = [matBaseName darkNum{ii} '.tiff'];
end

material.num        = 1;
material.latParm    = {5.411102};
material.structure  = 2;
material.modelPeaks = 1:10;
material.rList      = [ ...
    298 326; ...
    347 375; ...
    498 526; ...
    588 616; ...
    615 643; ...
    714 742; ...
    781 809; ...
    802 830; ...
    882 910; ...
    938 966; ...
    ].*instr.pixelSize;

cakeParms.bins(1)   = 36;%  ---> number of azimuthal bins
cakeParms.bins(2)   = 1800;%---> number of radial bins 1800
cakeParms.bins(3)   = 20;%  ---> number of angular bins 75
cakeParms.origin(1) = 1026.919; % ---> x center in pixels, fit2d convention
cakeParms.origin(2) = 1021.714; % ---> y center in pixels, fit2d convention
cakeParms.sector(1) = -360/cakeParms.bins(1)/2;% ---> start azimuth (min edge of bin) in degrees
cakeParms.sector(2) = 360-360/cakeParms.bins(1)/2;% ---> stop  azimuth (max edge of bin) in degrees
cakeParms.sector(3) = 250;  % ---> start radius (min edge of bin) in pixels
cakeParms.sector(4) = 1000; % ---> stop  radius (max edge of bin) in pixels
cakeParms.dims(1)   = instr.detectorSize/instr.pixelSize;%  ---> total number of rows in the full image
cakeParms.dims(2)   = instr.detectorSize/instr.pixelSize;%  ---> total
cakeParms.azim      = 0:360/cakeParms.bins(1):cakeParms.sector(2);
%%%%

%%%% apply geometric model
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
    darkData    = fliplr(darkData);
end

fileStps    = length(fileNum);
for i = 1:1:fileStps
    file.name   = fileList{i};
    [a b c]     = fileparts(file.name);
    instr.fileName  = [file.path b];
    
    if rebinImg
        [a b c]         = fileparts(file.name);
        instr.fileName  = [file.path b];
        
        imageData	= imread([instr.fileName '.tiff']);
        
        % FOR CHESS GE
        imdata      = imrotate(imageData, 180)';
        %%% FOR DISPLAY ONLY
        figure(101)
        imagesc(log(double(imdata-darkData))), axis square
        
        % FOR CHESS GE
        imageData   = fliplr(imageData);
        imageData   = imageData - darkData;
        
        % figure(101)
        % clf
        % imagesc(log(double(imageData)))
        % axis square tight
    end
    
    if useCeO2instr
        load([file.path CeO2List{i}])
        instr           = par.instr;
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
        return
        polImg  = reBin(mesh,instr,cakeParms,imageData);
        
        save([file.path b '.polImg.mat'],'polImg');
    else
        instr.fileName=[file.path b];
        disp([file.path b '.polImg.mat']);
        load([file.path b '.polImg.mat']);
    end
    
    figure(2)
    imagesc(log(polImg.intensity)), axis square
    
    if refitPeaks
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
            BuildESG(instr, cakeParms, polImgCor)
        end
        
        par.instr       = instr;
        par.material    = material;
        par.cakeParms   = cakeParms;
        par.peakInfo    = peakInfo;
        par.polImgCor   = polImgCor;
        par.qual        = qual;
        
        %%% AT THIS POINT, "par" VARIABLE IS SAVED
        %%% IN EACH FILE, HAVE CAKED IMAGES AND PEAK FITS
        save([instr.fileName '.par.mat'],'par');
    end
end
disp('caking and peak fitting complete ...')

%%% GENERATE SCATTERING VECTORS AND PLOT
for ihkl = 1:1:length(material.modelPeaks)
    q{ihkl}     = [];
    data{ihkl}  = [];
    for i = 1:1:fileStps
        file.name   = fileList{i};
        [a b c]     = fileparts(file.name);
        load([file.path b '.par.mat'])
        
        omega   = par.instr.omega;
        chi     = par.instr.chi;
        eta     = par.polImgCor.azimuth;
        
        Romega  = [ ...
            cosd(omega) 0 sind(omega); ...
            0 1 0; ...
            -sind(omega) 0 cosd(omega); ...
            ]';
        Rchi    = [ ...
            1 0 0; ...
            0 cosd(chi) -sind(chi); ...
            0 sind(chi) cosd(chi); ...
            ]';
        Rl2s    = Rchi*Romega;   %% LAB2SAM TRANSFORMATION
        
        %%% CONSTRUCT qj FOR EACH HKL / AZI
        tth	= atand(par.peakInfo.loc.*par.instr.pixelSize/par.instr.distance);
        th  = tth./2;
        
        data{ihkl}  = [data{ihkl}; par.peakInfo.amp(:,ihkl)];
        
        Ry  = [...
            cosd(th(ihkl)) 0 sind(th(ihkl)); ...
            0 1 0; ...
            -sind(th(ihkl)) 0 cosd(th(ihkl)); ...
            ];
        
        qj  = zeros(3, cakeParms.bins(1));
        q0  = [1 0 0]';
        for k = 1:1:cakeParms.bins(1)
            eta = par.cakeParms.azim(k);
            Rz  = [ ...
                cosd(eta) -sind(eta) 0; ...
                sind(eta) cosd(eta) 0; ...
                0 0 1; ...
                ];
            qj(:,k) = Rz*Ry*q0;
        end
        q{ihkl} = [q{ihkl} Rl2s*qj];
    end
    
    figure,
    PlotSPF(q{ihkl}', data{ihkl})
end