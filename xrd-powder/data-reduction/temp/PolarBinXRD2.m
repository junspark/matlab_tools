function polImg = PolarBinXRD2(mesh, instr, cakeParms, img)
Ifunc   = double(img);

L   = instr.detectorsize/instr.pixelsize;

% !!! THESE ARE IN THE CARTESIAN FRAME !!!
x0  = cakeParms.origin(1);   % in pixels
y0  = cakeParms.origin(2);   % in pixels

numEta  = cakeParms.bins(1);
numRho  = cakeParms.bins(2);

startRho    = cakeParms.sector(3);
stopRho     = cakeParms.sector(4);
dRho        = (stopRho - startRho)/numRho;

startEta    = cakeParms.sector(1);
stopEta     = cakeParms.sector(2);
deta        = (stopEta - startEta)/numEta;

R   = startRho:dRho:stopRho;
R   = repmat(R, cakeParms.bins(3)+1, 1)';

Rmean   = mean(R,2);
Rlist   = (Rmean(1:end-1) + Rmean(2:end))./2;

polImg.radius    = zeros(size(Rmean,1)-1, numEta);
polImg.azimuth   = cakeParms.azim;
polImg.intensity = zeros(size(Rmean,1)-1, numEta);
for ii=1:numEta
    msg = sprintf('Processing sector %d of %d', [ii, numEta]);
    fprintf(1, 'INFO: %s\n', msg);
    tic;
    
    eta         = polImg.azimuth(ii);
    eta_ini     = eta - deta/2;
    eta_fin     = eta + deta/2;
    eta_step    = cakeParms.bins(3);
    eta_grid    = linspace(eta_ini, eta_fin, eta_step + 1);
    TH          = repmat(eta_grid, numRho+1, 1);
    
    [x, y]  = pol2cart(deg2rad(TH), R);
    
    X   = x0 + x;
    Y   = y0 + y;
    
    figure(1)
    plot(X, Y,'b.')
    V   = zeros(size(X));
    for i=1:size(X,1)
        for j=1:size(X,2)
            XY=[X(i,j);Y(i,j)];
            plot(XY(1), XY(2), 'ko')
            V(i,j) = DataCoordinates(XY, L, mesh, Ifunc);
        end
    end
    % Ilist   = BuildMeshPolarXRD(R, V, 360/cakeParms.bins(1)/numEta, mesh.qrule);
    Ilist   = BuildMeshPolarXRD(R, V, eta_grid(2) - eta_grid(1), mesh.qrule);
    keyboard
    figure,
    plot(Rlist, Ilist, '.')
    polImg.radius(:,ii)    = Rlist';
    polImg.intensity(:,ii) = Ilist';
    time20(ii)=toc;
    disp(sprintf('Processing time for sector %d is %1.4f', [ii, time20(ii)]));
end

polImg.radius=polImg.radius';
polImg.intensity=polImg.intensity';