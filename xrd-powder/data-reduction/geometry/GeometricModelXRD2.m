function mappedData = GeometricModelXRD2(centTrans, D, ...
    gammaYprime, gammaXhatPrime, ...
    data, azim, DistPars)

% GeometricModelXRD2 - Apply geometric model and detector distortion to the
% raw 2theta fits
%
% USAGE:
%
% INPUTS:
%
% 1) centTrans
% 2) D
% 3) gammaYprime
% 4) gammaXhatPrime
% 5) data
% 6) azim
% 7) DistType
% 8) DistPars

numSpectra = size(data,2);
if length(azim) ~= numSpectra
    error('length of azimuth list does not correspond to the # of spectra');
end

%% origins of the scattering center and detector frames
P1 = [0, 0, 0]';
P2 = [0, 0, D]';

%% basis vectors for the tilt-corrected (primed) frame
Xprime = [1 0 0]';
Yprime = [0 1 0]';
Zprime = [0 0 1]';

%% rotation 1: gammaYprime about Yprime (fit2d ROTX?)
RYprime = RMatOfExpMap(gammaYprime * Yprime);

%% the transformed Xprime axis...
XhatPrime = RYprime * Xprime;

%% rotation 2: gammaXhatPrime about XhatPrime (fit2d ROTY?)
RXhatPrime = RMatOfExpMap(gammaXhatPrime * XhatPrime);

%% change of basis matrix from hatPrime to Prime
RHatPrimeToPrime = RXhatPrime * RYprime;

%% assign origin displacement
tHatPrimeToHatdPrime = centTrans;
tUnToPrime = P2;
mappedData=[];
for i = 1:numSpectra
    numPts = size(data(:,i), 1);
    
    %% assign radial coordinates and associated intensities
    %% PROPOSED BY Lee J.H.
    f = DistPars(1)*((data(:,i)/DistPars(5)).^(DistPars(3)))*cosd(4*azim(i))...
        + DistPars(2)*((data(:,i)/DistPars(5)).^DistPars(4)) + 1;
    rhoHatdPrime  = f.*data(:,i);
    
    %% assign azimuth (note units!)
    etaHatdPrime = deg2rad(azim(i));
    
    %% Convert to cartesian coord's in integration (hat double-prime) frame
    xHatdPrime = rhoHatdPrime * cos(etaHatdPrime);
    yHatdPrime = rhoHatdPrime * sin(etaHatdPrime);
    
    %% apply center translation on to get components in IP-centered (hat
    %% prime) frame
    xHatPrime = xHatdPrime + repmat(tHatPrimeToHatdPrime(1), [numPts, 1]);
    yHatPrime = yHatdPrime + repmat(tHatPrimeToHatdPrime(2), [numPts, 1]);
    
    rhoHatPrime = sqrt(xHatPrime.^2 + yHatPrime.^2);
    etaHatPrime = rad2deg(atan2(yHatPrime, xHatPrime));
    
    %% full components in IP-centered (hat prime) frame
    P3HatPrime = [xHatPrime, yHatPrime, zeros(numPts, 1)]';
    
    %% rotate components into tilt-corrected (prime) frame
    P3Prime    = RHatPrimeToPrime * P3HatPrime;
    
    %% apply translation to get to sample (unprimed) frame
    P3 = P3Prime + repmat(tUnToPrime, [1, numPts]);
    
    %% get two-theta from dot products
    dotProds = dot(...
        UnitVector(P3 - repmat(P1, [1, numPts])), ...
        repmat(UnitVector(P2 - P1), [1, numPts]), 1);
    
    %% two-theta
    measTTH = acos(dotProds);
    
    %% transform data
    tmpData     = rad2deg(measTTH(:));
    mappedData  = [mappedData tmpData];
end