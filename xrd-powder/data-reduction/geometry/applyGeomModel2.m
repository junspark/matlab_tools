%%% LOAD DARK
if rebinImg
    fileStps    = length(darkNum);
    for i = 1:1:fileStps
        filename    = darkList{i};
        [a b c]     = fileparts(filename);
        fileName    = [file.path b];
        if i == 1
            darkData    = imread([fileName '.tiff']);
        else
            darkData    = darkData + imread([fileName '.tiff']);
        end
    end
    darkData    = darkData./fileStps;
    
    % FOR CHESS GE
    darkData    = fliplr(darkData);
end

%%% LOAD INSTR PRM
energy      = 0;
distance	= 0;
centers     = [0 0];
gammaX      = 0;
gammaY      = 0;
radOffset   = 0;

DistType    = 0;

a0	= 0;
a1  = 0;
a2  = 0;
n1  = 0;
n2  = 0;
t1  = 0;
e1  = 0;

for i = 1:1:length(CeO2fileNum)
    filename    = CeO2List{i};
    [a b c]     = fileparts(filename);
    fileName    = [file.path b, '.mat'];
    load(fileName)
    
    energy      = par.instr.energy + energy;
    distance    = par.instr.distance + distance;
    centers     = par.instr.centers + centers;
    gammaX      = par.instr.gammaX + gammaX;
    gammaY      = par.instr.gammaY + gammaY;
    radOffset   = par.instr.radOffset + radOffset;
    
    DistType    = par.instr.DistType + DistType;
    
    a0  = par.instr.a0 + a0;
    a1  = par.instr.a1 + a1;
    a2  = par.instr.a2 + a2;
    n1  = par.instr.n1 + n1;
    n2  = par.instr.n2 + n2;
    t1  = par.instr.t1 + t1;
    e1  = par.instr.e1 + e1;
end

instr.energy        = energy./length(CeO2fileNum);
instr.wavelength    = keV2Angstrom(instr.energy);% wavelength (Angstrom)
instr.distance      = distance./length(CeO2fileNum);
instr.centers       = centers./length(CeO2fileNum);
instr.gammaX        = gammaX./length(CeO2fileNum);
instr.gammaY        = gammaY./length(CeO2fileNum);
instr.radOffset     = radOffset./length(CeO2fileNum);
instr.rock          = 5;
instr.tthRng        = [3 14];
instr.StepSize      = .5*100;
instr.WindowSize    = 100;
instr.Base          = 2;
instr.Levels        = 12;
instr.FWHHFilter    = .04;
instr.HeightFilter  = 900;
instr.PeakLocation  = .5;
instr.DistType      = DistType;

instr.a0    = a0./length(CeO2fileNum);
instr.a1    = a1./length(CeO2fileNum);
instr.a2    = a2./length(CeO2fileNum);
instr.n1    = n1./length(CeO2fileNum);
instr.n2    = n2./length(CeO2fileNum);
instr.t1    = t1./length(CeO2fileNum);
instr.e1    = e1./length(CeO2fileNum);

for i = 1:1:length(matfileNum)
    if rebinImg
        [a b c]     = fileparts(fileList{i});
        fileName    = [file.path b];
        
        instr.fileName  = fileName;
        
        imageData	= imread([fileName '.tiff']);
        
        % FOR CHESS GE
        imageData   = fliplr(imageData);
        
        imageData   = imageData - darkData;
        
        if i == 1
            mesh    = buildMeshXRD(imageData);
        end
        polImg      = reBin(mesh, instr, cakeParms, imageData);
    end
    
    instr.omega = stndOmega(i);
    instr.chi   = stndChi(i);
    
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
    j   = length(DistPars);
    
    x(1)        = instr.centers(1);
    x(2)        = instr.centers(2);
    x(3)        = instr.gammaY;
    x(4)        = instr.gammaX;
    x(5)        = instr.distance;
    x(6:5+j)    = DistPars;
    
    %%% CORRECT BINNED IMAGE FOR INSTR GEOMETRY
    DATA    = cell(1,cakeParms.bins(1));
    for ii = 1:cakeParms.bins(1)
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
    par.polImgCor   = polImgCor;
    
    save([fileName '.par.mat'],'par');
end