function polImg = PolarBinXRD(mesh, instr, cakeParms, img)
Ifunc   = double(img);

L   = instr.detectorsize/instr.pixelsize;

% !!! THESE ARE IN THE CARTESIAN FRAME !!!
x0  = cakeParms.origin(1);   % in pixels
y0  = cakeParms.origin(2);   % in pixels

startRho    = cakeParms.sector(3);
stopRho     = cakeParms.sector(4);

numEta = cakeParms.bins(3);
numRho = cakeParms.bins(2);

R   = repmat(cakeParms.sector(3):(cakeParms.sector(4)-cakeParms.sector(3))/numRho:cakeParms.sector(4),numEta+1,1)';

Rlist=zeros(size(R,2)-1,1);
for jj=1:size(R,1)-1
    Rval=zeros(size(R,2)-1,1);
    for ii=1:size(R,2)-1
        Rval(ii)=mean([R(jj,ii) R(jj+1,ii)]);
    end
    
    Rlist(jj)=mean(Rval);
end

polImg.radius    = zeros(size(R,1)-1,cakeParms.bins(1));
polImg.azimuth   = cakeParms.azim;
polImg.intensity = zeros(size(R,1)-1,cakeParms.bins(1));

for ii=1:cakeParms.bins(1)
    msg = sprintf('Processing sector %d of %d', [ii, cakeParms.bins(1)]);
    fprintf(1, 'INFO: %s\n', msg);
    tic;
    TH=repmat(polImg.azimuth(ii)-360/cakeParms.bins(1)/2:360/cakeParms.bins(1)/numEta:polImg.azimuth(ii)+360/cakeParms.bins(1)/2,size(R,1),1);
    [x,y]=pol2cart(deg2rad(TH),R);
    X=x0+x;
    Y=y0+y;
  
    V=zeros(size(X));
    for i=1:size(X,1)
        for j=1:size(X,2)
            XY=[X(i,j);Y(i,j)];
            V(i,j) = DataCoordinates(XY,L, mesh, Ifunc);
        end
    end
    
    Ilist=buildMeshPolarXRD(R,V,360/cakeParms.bins(1)/numEta,mesh.qrule);
    
    polImg.radius(:,ii)    = Rlist';
    polImg.intensity(:,ii) = Ilist';
    time20(ii)=toc;
    msg2 = sprintf('Processing time for sector %d is %1.4f', [ii, time20(ii)]);
end

polImg.radius=polImg.radius';
polImg.intensity=polImg.intensity';