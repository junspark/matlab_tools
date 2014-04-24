function polImg = PolarBinXRD2(mesh, instr, cakeParms, img)
Ifunc   = double(img);

L   = instr.detectorsize/instr.pixelsize;

% !!! THESE ARE IN THE CARTESIAN FRAME !!!
x0  = cakeParms.origin(1);   % in pixels
y0  = cakeParms.origin(2);   % in pixels

startAzi    = cakeParms.sector(1);
stopAzi     = cakeParms.sector(2);
numAzi      = cakeParms.bins(1);
dAzi        = (stopAzi - startAzi) / numAzi;

startRho    = cakeParms.sector(3);
stopRho     = cakeParms.sector(4);
numRho      = cakeParms.bins(2);

numEta = cakeParms.bins(3);

R       = linspace(startRho, stopRho, numRho + 1);
R       = repmat(R, numEta+1, 1)';
Rlist   = (sum(R(1:end-1,:),2) + sum(R(2:end,:),2))./(numEta + 1)/2;

polImg.radius    = zeros(size(R,1)-1,cakeParms.bins(1));
polImg.azimuth   = cakeParms.azim;
polImg.intensity = zeros(size(R,1)-1,cakeParms.bins(1));
for ii=1:cakeParms.bins(1)
    msg = sprintf('Processing sector %d of %d', [ii, cakeParms.bins(1)]);
    fprintf(1, 'INFO: %s\n', msg);
    tic;
    
    eta         = polImg.azimuth(ii);
    eta_ini     = eta - dAzi/2;
    eta_fin     = eta + dAzi/2;
    eta_grid    = linspace(eta_ini, eta_fin, numEta + 1);
    eta_grid    = repmat(eta_grid, numRho + 1, 1);

    [x, y]  = pol2cart(deg2rad(eta_grid),R);
    
    X   = x0 + x;
    Y   = y0 + y;
    
    figure(1)
    plot(X, Y, 'k:')
    
    V=zeros(size(X));
    for i=1:size(X,1)
        for j=1:size(X,2)
            XY=[X(i,j);Y(i,j)];
            V(i,j) = DataCoordinates(XY, L, mesh, Ifunc);
        end
    end
    % keyboard
    Ilist   = BuildMeshPolarXRD(R, V, 360/cakeParms.bins(1)/numEta, mesh.qrule);
    polImg.radius(:,ii)    = Rlist';
    polImg.intensity(:,ii) = Ilist';
    timestamp   = toc;
    disp(sprintf('Processing time for sector %d is %1.4f', [ii, timestamp]));
end

polImg.radius=polImg.radius';
polImg.intensity=polImg.intensity';