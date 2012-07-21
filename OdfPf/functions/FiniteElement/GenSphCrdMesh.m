function smesh = GenSphCrdMesh(ntheta, nphi, thetaStart, thetaEnd, phiStart, phiEnd)
% GENSPHCRDMESH - Generate a hemisphere mesh based on spherical coordinates.
%
%   smesh = GenSphCrdMesh(ntheta, nphi, thetaStart, thetaEnd, phiStart, phiEnd)
%
%   ntheta     is 1 x 1, the number of subdivisions in azimuth
%   nphi       is 1 x 1, the number of subdivisions in declination
%   thetaStart is 1 x 1, the starting value of theta (0, 2*pi)
%   thetaEnd   is 1 x 1, the final value of phi (0, 2*pi) > thetaStart
%   phiStart   is 1 x 1, the starting value of phi (0, pi)
%   phiEnd     is 1 x 1, the final value of phi (0, pi) > phiStart
%
%   smesh is a MeshStructure on the hemisphere.
%
%   Notes:
%
%   *)  This has been modified to pass ranges in theta and phi.
%
%   *)  It could also be modified to generate equivalences for
%       the projective sphere.
%
tdiv = thetaStart:(thetaEnd - thetaStart)/ntheta:thetaEnd;
pdiv = phiStart:(phiEnd - phiStart)/nphi:phiEnd;
%
[phi, theta] = meshgrid(pdiv, tdiv);

wtht    = repmat(0.05*phi(1, :), [ntheta + 1, 1]) + theta;
fullCon = delaunay(phi(:), wtht(:));
nel     = size(fullCon, 1);

crds = [phi(:), theta(:)];

overlapChk = sum(ismember([0, 2*pi], tdiv)) == 2;

newCon = fullCon;

if overlapChk
  overlap360 = (ntheta + 1)*[1:nphi+1];
  mapTo      = overlap360 - ntheta;
  for i = 1:length(mapTo)
    fullCon(find(fullCon == overlap360(i))) = mapTo(i);
    %crds = overlap360(i);
  end

  indexi = cell(1, nphi);
  for i = 1:nphi
    indexi{i} = find(fullCon >= i*(ntheta + 1) & fullCon <= (i + 1)*(ntheta + 1));
  end

  %% newCon = fullCon;
  for i = 1:size(indexi, 2)
    newCon(indexi{i}) = fullCon(indexi{i}) - i;
  end

  newphi   = phi(1:end - 1, :);
  newtheta = theta(1:end - 1, :);

  newcrd = [newphi(:), newtheta(:)];
end

if phiStart == 0 & overlapChk  
  indxr = [1:ntheta];

  subtractMe = ~ismember(newCon, indxr);
  degenChk = sum(~subtractMe, 2) == 2;

  newCon(~subtractMe) = 1;
  newCon(subtractMe)  = newCon(subtractMe) - indxr(end) + 1;
  newCon = newCon(~degenChk, :);

  newphi   = newphi(indxr(end):end);
  newtheta = [newtheta(1), newtheta(indxr(end) + 1:end)];

  newcrd = [newphi(:), newtheta(:)];
elseif phiStart == 0 & ~overlapChk
  indxr = [1:ntheta + 1];

  subtractMe = ~ismember(fullCon, indxr);
  degenChk = sum(~subtractMe, 2) == 2;

  %% newCon = fullCon;
  newCon(~subtractMe) = 1;
  newCon(subtractMe)  = fullCon(subtractMe) - indxr(end) + 1;
  newCon = newCon(~degenChk, :);

  newphi   = phi(indxr(end):end);
  newtheta = [theta(1), theta(indxr(end) + 1:end)];

  newcrd = [newphi(:), newtheta(:)];
elseif phiStart ~= 0 & ~overlapChk
  newcrd = crds;
  %% newCon = fullCon;
end

%% IF MESH CONTAINS NEGATIVE POLE
if phiEnd == pi
  maxpt = max(max(newCon));
  indxr = [maxpt-ntheta:maxpt];

  subtractMe = ~ismember(newCon, indxr);
  degenChk = sum(~subtractMe, 2) == 2;

  newCon(~subtractMe) = indxr(1);
  newCon = newCon(~degenChk, :);

  newphi   = newphi(1:indxr(1));
  newtheta = newtheta(1:indxr(1));

  newcrd = [newphi(:), newtheta(:)];
end

[X, Y, Z] = sph2cart(newcrd(:, 2), pi/2 - newcrd(:, 1), 1);
newcrd_c = [X, Y, Z]';

smesh = MeshStructure(newcrd_c, newCon');
