function mappedData = XFormRadialSpectra(centTrans, D, gammaYprime, gammaXhatPrime, data, azim, tthRange,DistPars,DistType)
% XFORMRADIALSPECTRA - Generate intensity vs. 2-theta spectra by applying
% geometric transform to 1-d radial spectra from area detector images.
%
% USAGE:
%
% mappedData = XFormRadialSpectra(...
% 	centTrans, D, ...
% 	gammaYprime, gammaXhatPrime, ...
% 	radOffset, ...
% 	data, azim, ...
% 	tthMin, tthMax)
%
% INPUTS:
%
% 1) centTrans
% 2) D
% 3) gammaYprime
% 4) gammaXhatPrime
% 5) radOffset
% 6) data
% 7) azim
% 8) tthRange

numSpectra = length(data);
if length(azim) ~= numSpectra
    error('length of azimuth list does not correspond to the # of spectra');
end

%% origins of the scattering center and detector frames
P1 = [0, 0, 0]';
P2 = [0, 0, -D]';

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

for i = 1:numSpectra
    numPts = size(data{i}, 1);
    
    %% assign radial coordinates and associated intensities; also apply
    %% radial offset (e.g. for correction of mar345 data)
    if DistType==0
        rhoHatdPrime  = data{i}(:, 1);
    elseif DistType==1
        rhoHatdPrime  = data{i}(:, 1) + DistPars(1);
    elseif DistType==2
        f = DistPars(1)*((data{i}(:, 1)/DistPars(7)).^(DistPars(3)))*cosd(DistPars(5)*azim(i)+DistPars(6))...
            + DistPars(2)*((data{i}(:, 1)/DistPars(7)).^DistPars(4)) + 1;
        rhoHatdPrime=f.*data{i}(:, 1);
    elseif DistType==3
        f = DistPars(1)*((data{i}(:, 1)/DistPars(7)).^(DistPars(3)))*cosd(DistPars(5)*azim(i)+DistPars(6))...
            + DistPars(2)*((data{i}(:, 1)/DistPars(7)).^DistPars(4)) + 1;
        rhoHatdPrime=f.*data{i}(:, 1);
    end
    
    measIntensity = data{i}(:, 2);
    
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
    tmpData = [rad2deg(measTTH(:)), measIntensity(:), etaHatPrime(:)];

    if nargin > 7
        tthRange = sort(tthRange);
        tthMin = tthRange(1);
        tthMax = tthRange(2);

        tmpData = tmpData((tmpData(:, 1) >= tthMin & tmpData(:, 1) <= tthMax), :);
    end

    mappedData{i} = tmpData;
end
