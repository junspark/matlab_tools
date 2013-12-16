clear all
close all
clc

%%% DATA REDUCTION FLAGS
UseCeO2InstPrms = 1;    % USE DETECTOR PARAMETERS OBTAINED FROM CALIBRANT
% MakePNGs        = 0;    % OUTPUT PNG FILES OF IMAGES / FRAMES
BinImage        = 1;    % BIN IMAGE
FitPeaks        = 1;    % FIT PEAKS
GenerateESG     = 1;    % GENERATES ESG FOR TEXTURE ANALYSIS
DetOrientation  = 2;    % GE NUMBER (0 - 4 APS GE)
OutputPath      = 'C:\Users\parkjs\Documents\Projects\TiNi\schaffer_201306';

%%% DEFINE IMAGE SERIES
%%% IF CALIBRANT IMAGES, SET VARIABLE UseCeO2InstPrms TO 0 TO OBTAIN
%%% CALIBRATION AND SAVE TO DETPAR FILE
%%% IF MATERIAL DATA, SET VARIABLE UseCeO2InstPrms TO 1 AND PROGRAM WILL
%%% LOOK FOR DETPAR FILE
IMSeries.pname  = 'W:\schaffer_201306';
IMSeries.froot  = 'CeO2_';    % NiTi_sam23_ FOR MATERIAL, CeO2_FOR CEO2
IMSeries.fnumi  = 679;       %20:1:638 FOR MATERIAL, 679 FOR CEO2
IMSeries.fnumf  = 679;
IMSeries.finc   = 1;
IMSeries.ndigit = 5;
IMSeries.fnum   = IMSeries.fnumi:IMSeries.finc:IMSeries.fnumf;
IMSeries.fext   = 'ge2';
for i = 1:1:length(IMSeries.fnum)
    IMSeries.flist{i,1} = [IMSeries.froot, ...
        sprintf(['%0', num2str(IMSeries.ndigit), 'd'], IMSeries.fnum(i)), ...
        '.', IMSeries.fext];
end
IMSeries.detpars    = 'detector_hexrd.detpar';

IMSeries.dark_pname = 'W:\schaffer_201306';
IMSeries.dark_froot = 'dark_';
IMSeries.dark_fnum  = [680; 680];
for i = 1:1:length(IMSeries.dark_fnum)
    IMSeries.dark_flist{i,1}    = [IMSeries.dark_froot, ...
        sprintf(['%0', num2str(IMSeries.ndigit), 'd'], IMSeries.dark_fnum(i)), ...
        '.', IMSeries.fext];
end

IMSeries.omega  = 0;

%%% INSTRUMENT PARAMSTERS
Instr.energy        = 86;                           % keV
Instr.wavelength    = keV2Angstrom(Instr.energy);   % wavelength (Angstrom)
Instr.distance      = 1490;                         % mm
Instr.cenX          = 204.8;            % mm
Instr.cenY          = 204.8;            % mm
Instr.gammaX        = 0;                % deg  
Instr.gammaY        = 0;                % deg
Instr.detX          = 409.6;    % mm
Instr.detY          = 409.6;    % mm
Instr.pixelSizeX    = 0.200;    % mm
Instr.pixelSizeY    = 0.200;    % mm
Instr.rock          = 0;

%%% DETECTOR DISTORTION PARAMETERS
Instr.p0	= 0;
Instr.p1    = 0;
Instr.p2    = 0;
Instr.n0    = 2;
Instr.n1    = 2;
Instr.n2    = 2;

%%% MATERIAL PARAMETERS
material(1).LattParm    = 5.411102;
material(1).Structure   = 2;
material(1).ModelPeaks  = 1:1:10;

%%% CAKING PARAMETERS
cakeParms.bins(1)   = 36;       % ---> number of azimuthal bins
cakeParms.bins(2)   = 1800;     % ---> number of radial bins 1800
cakeParms.bins(3)   = 75;       % ---> number of angular bins 75
cakeParms.origin(1) = 1026.919; % ---> x center in pixels
cakeParms.origin(2) = 1021.714; % ---> y center in pixels
cakeParms.sector(1) = -360/cakeParms.bins(1)/2;     % ---> start azimuth (min edge of bin) in degrees
cakeParms.sector(2) = 360-360/cakeParms.bins(1)/2;  % ---> stop  azimuth (max edge of bin) in degrees
cakeParms.sector(3) = 250;  % ---> start radius (min edge of bin) in pixels
cakeParms.sector(4) = 1000; % ---> stop  radius (max edge of bin) in pixels
cakeParms.dims(1)   = Instr.detX/Instr.pixelSizeX;   %  ---> total number of rows in the full image
cakeParms.dims(2)   = Instr.detY/Instr.pixelSizeY;   %  ---> total
cakeParms.azim      = 0:360/cakeParms.bins(1):cakeParms.sector(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if UseCeO2InstPrms
    pfname  = fullfile(OutputPath, IMSeries.detpars);
    detpars = load(pfname);
    
    Instr.cenX      = detpars(1,1);
    Instr.cenY      = detpars(2,1);
    Instr.distance  = detpars(3,1);
    Instr.gammaX	= detpars(4,1);
    Instr.gammaY    = detpars(5,1);
    
    Instr.p0	= detpars(6,1);
    Instr.p1    = detpars(7,1);
    Instr.p2    = detpars(8,1);
    Instr.n0    = detpars(9,1);
    Instr.n1    = detpars(10,1);
    Instr.n2    = detpars(11,1);
end

if BinImage
    %%% GENERATE DARK FILE FIRST
    ct  = 0;
    imdata_dark = zeros(cakeParms.dims(1), cakeParms.dims(2));
    for i = 1:1:length(IMSeries.dark_fnum)
        pfname  = fullfile(IMSeries.pname, IMSeries.dark_flist{i,1});
        
        d       = dir(pfname);
        nframes = (d.bytes - 8192)/2048/2048/2;
        for j = 1:1:nframes
            imdata_dark = imdata_dark + double(NreadGE(pfname,j));
            ct  = ct + 1;
        end
    end
    imdata_dark = imdata_dark./ct;
    
    figure(1000)
    imagesc(imdata_dark)
    title('dark field - may not be in correct orientation')
    axis equal tight
    hold off
    
    %%% LOOK AT INDIVIDUAL IMAGE STACK 
    for i = 1:1:length(IMSeries.fnum)
        pfname  = fullfile(IMSeries.pname, IMSeries.flist{i,1});
        
        d       = dir(pfname);
        nframes = (d.bytes - 8192)/2048/2048/2;
        for j = 1:1:nframes
            imdata	= double(NreadGE(pfname,1)) - imdata_dark;
            idx     = imdata < 0;
            imdata(idx) = 0;
            
            %%% FOR VIEWING
            %%% H IS XL, V IS YL
            imdata  = rot90(imdata, 1);
            
            figure(1000)
            % imagesc(log(abs(imdata)))
            imagesc(imdata, [0 1500])
            title([IMSeries.flist{i,1}, ' - image data - ', num2str(i)], 'Interpreter', 'none')
            axis equal tight
            colorbar vert
            hold off
            
            return
            pause(0.05)
        end
    end
end
return
%%%%

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