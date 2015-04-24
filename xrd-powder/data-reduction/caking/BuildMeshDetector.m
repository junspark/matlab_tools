function mesh = BuildMeshDetector(Lx, cakeParms)
% BuildMeshDetector - Build mesh crd and con for an XRD image
%   
%   mesh = BuildMeshDetector(Lx, cakeParms)
%
%   Input:
%
%   Lx
%       number of pixels in the image in the horizontal direction
%
%   cakeParms
%       caking parameters 
%
%   Output:
%   
%   mesh 
%       mesh structure of the detector with appropriate fields
%
%   Note: 
%       Only works with square images for now (2015-04-24)


%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD STUFF
% function mesh = BuildMeshDetector(Lx, Ly, cakeParms)
%   Ly 
%       number of pixels in the image in the vertical direction
%


%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE MESH
%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin == 1
%     disp('square image')
%     Ly = Lx;
% elseif nargin == 2
%     disp('rectangular image')
% end
Ly = Lx;

numel   = (Lx - 1) * (Ly - 1) * 2;
con     = zeros(3, numel);

for i = 1:1:(Ly - 1)
    idx1    = 1 + 2*(Lx - 1) * ( i - 1 );
    idx2    = 2*(Lx - 1) * i;
    
    xL  = Lx*(i - 1)+1;
    xR  = Lx*i;
    xLR = (xL + 1):1:(xR - 1);
    xLR = repmat(xLR, 2, 1);
    
    con(1, idx1: idx2)  = [xL xLR(:)' xR];
    
    xL  = Lx*(i - 1) + 2 : 1 : Lx * i;
    xR  = Lx*i + 2 : 1 : Lx * (i + 1);
    xLR = [xL; xR];
    
    con(2, idx1: idx2)  = xLR(:)';
    
    xL  = Lx * i + 1;
    xR  = Lx * (i + 1) - 1;
    xLR = xL:xR;
    xLR = repmat(xLR, 2, 1);
    con(3, idx1: idx2)  = xLR(:)';
    
end
x   = repmat(1:Lx, 1, Ly);
y   = repmat(1:Ly, Lx, 1);
y   = y(:)';

crd = [x; y];

load('qr_trid03p06.mat');
mesh    = MeshStructureXRD(crd, con, numel, qr_trid03p06);

%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE LOOKUP TABLE FOR BINNING
L   = fix(Lx);  % ASSUME SQUARE IMAGE FOR NOW

% !!! THESE ARE IN THE CARTESIAN FRAME !!!
x0  = cakeParms.origin(1);   % in pixels
y0  = cakeParms.origin(2);   % in pixels

startAzi    = cakeParms.sector(1);
endAzi      = cakeParms.sector(2);

startRho    = cakeParms.sector(3);
endRho      = cakeParms.sector(4);

%%% NUMBER OF AZIMUTHAL BINS OVER ANGULAR RANGE DEFINED BY cakeParms.sector(1) AND cakeParms.sector(2)
numAzi  = cakeParms.bins(1);
%%% NUMBER OF RADIAL POINTS PER AZIMHUTHAL BIN OVER RADIAL RANGE DEFINED BY startRho AND endRho
numRho  = cakeParms.bins(2);
%%% NUMBER OF ETA POINTS PER AZIMUTH
numEta  = cakeParms.bins(3);

dAzi    = (endAzi - startAzi)/numAzi;
dEta    = dAzi/numEta;

R   = startRho:(endRho-startRho)/numRho:endRho;
R   = repmat(R, numEta + 1, 1)';

for ii = 1:1:numAzi
    fprintf('Producing lookup table for sector %d of %d\n', ii, numAzi);
        
    azi_ini = cakeParms.azim(ii) - dAzi/2;
    azi_fin = cakeParms.azim(ii) + dAzi/2;
    
    TH  = azi_ini:dEta:azi_fin;
    TH  = repmat(TH, numRho + 1, 1);
    
    [x, y]	= pol2cart(deg2rad(TH),R);
    x   = x0 + x; 
    y   = y0 + y;
    
    tic;
    for i = 1:1:(numRho + 1)
        for j = 1:1:(numEta + 1)
            xy  = [x(i,j); y(i,j)];
            
            % square grid
            gridNum     = (L-1)*(fix(xy(2))-1)+fix(xy(1));
            elemNums    = [gridNum*2-1 gridNum*2];
            
            P   = crd(:, con(:,elemNums(1)));
            A1  = [P(1,1)-P(1,3) P(1,2)-P(1,3);P(2,1)-P(2,3) P(2,2)-P(2,3)];
            B1  = [xy(1)-P(1,3); xy(2)-P(2,3)];
            Y1  = A1\B1;
            
            if all(Y1>=0)
                fele    = elemNums(1);
                fcrd{ii,i,j}    = [Y1(1) Y1(2) 1-Y1(1)-Y1(2)];
            else
                P   = crd(:,con(:,elemNums(2)));
                A2  = [P(1,1)-P(1,3) P(1,2)-P(1,3);P(2,1)-P(2,3) P(2,2)-P(2,3)];
                B2  = [xy(1)-P(1,3); xy(2)-P(2,3)];
                Y2  = A2\B2;
                
                fele            = elemNums(2);
                fcrd{ii,i,j}    = [Y2(1) Y2(2) 1-Y2(1)-Y2(2)];
            end
            fcon{ii,i,j}    = con(:, fele);
        end
    end
    
    dtime   = toc;
    fprintf('Processing time for sector %d is %1.4f\n', ii, dtime);
end
mesh.fcon   = fcon;
mesh.fcrd   = fcrd;
