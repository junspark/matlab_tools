function q = MinRepOfEqvClass(qeqv, leftSymTag)
% MINREPOFEQVCLASS - Put equivalence class of quaternions under
%  left and right symmetries into fundamental region (disorientation space)
%
%   q = MisorientationOfEqvClass(qeqv, leftSymTag)
%
%   qeqv is 4 x l*m x n, an equivalence class for an array quaternions
%                        under left and right symmetries of order l and m,
%                        respectively with l <= m.
%   leftSymTag is 1 x 1, a string containing the name of the left symmetry
%                        group (see note)
%
%   q is 4 x n, the array of quaternions lying in the composite
%               fundamental region for the symmetry groups
%               in question
%
%   Notes:
%
%   11/19/2003 *)  This routine is typically called by Misorientation.m and
%                  ToFundamentalRegionQ.m, where symmetries have already
%                  been applied to the points.  It is intended to simply
%                  store the axial vector restriction conditions for
%                  obtaining the representative set of quaternions.  The
%                  right symmetry is as of yet not needed, as this
%                  typically is the group of greater or equal order and
%                  only restricts the maximum misorientation angle (dist
%                  from origin).
%
% The following follows the treatment of Heinz and Neumann (Acta Cryst.
% (1991). A47, 780-789) for selecting a representative (mis)orientation
% R under the influence of two lattice-lattice/sample-lattice symmetry
% groups, G1 and G2. The indicies i and j refer to the order of the
% symmetry groups G1 and G2, repsectively.  It should be noted that we
% adopt the convention in which as orientation describes the components
% of a crystal quantity vc in the sample reference (vs), hence in the
% presence of lattice symmetry G of order k we have
%
%   vs = R * (G * vc).
%
% Given a G1 and G2 where i = 1:n and j = 1:m, we obtain n x m
% equivalent descriptions of R:
%
% R_ij = G2_j * R * G1_i
%      = G2_j * R * G1_i * (G2_j * G2_j^-1)
%      = G2_j * R * (G1_i * G2_j) * G2_j^-1
%      = G2_j * (R * G1_i * G2_j) * G2_j^-1
%      = G2_j * Rmin * G2_j^-1
%
npts  = size(qeqv, 3); % number of points
neqv  = size(qeqv, 2); % size of equivalence class
reverseCols = neqv:-1:1;

% zero tolerance (could make this another optional input)
ztol = 1e-10;

% Fix the near-zero entries for comparison
qeqv(find(abs(qeqv) <= ztol)) = 0;

% find the points whose axial vector lies in the sterographic traingle
% of G2 based on leftSymTag
newq = zeros(4, npts);
for i = 1:npts
    %% Process cases for currently implemented symmetry groups
    if strcmp(leftSymTag, 'triclinic')
        % Admittedly, a little silly, but for the sake of consistency
        % and readability, I've let triclinic symmetry in the loop
        [qmax, imax] = max(abs(qeqv(1, :, i)));
        %
        newq(:, i) = qeqv(:, imax, i);
    elseif strcmp(leftSymTag, 'orthorhombic')
        % Using SST of [100], [010], [001]
        [qsort, isort] = sort(abs(qeqv(1, :, i)));
        tmp = qeqv(:, isort, i);
        tmp = tmp(:, reverseCols);
        %
        firstQuadrant  = find((tmp(2, :) >= 0) & (tmp(3, :) >= 0) & (tmp(4, :) >= 0));
        secondQuadrant = find((tmp(2, :) >= 0) & (tmp(3, :) <= 0) & (tmp(4, :) >= 0));
        %
        newq(:, i) = tmp(:, min(union(firstQuadrant, secondQuadrant)));
        %
    elseif strcmp(leftSymTag, 'tetragonal')
        % Using SST of [100], [110], [001]
        [qsort, isort] = sort(abs(qeqv(1, :, i)));
        tmp = qeqv(:, isort, i);
        tmp = tmp(:, reverseCols);
        tmp1 = single(tmp);
        atmp1 = single(abs(tmp));
        %
        firstQuadrant  = find((tmp1(2, :) >= tmp1(3, :)) & (tmp1(3, :) >= 0) & (tmp1(4, :) >= 0));
        secondQuadrant = find((tmp1(2, :) >= -tmp1(3, :)) & (tmp1(3, :) <= 0) & (tmp1(4, :) >= 0));
        %
        newq(:, i) = tmp(:, min(union(firstQuadrant, secondQuadrant)));
        %
    elseif strcmp(leftSymTag, 'cubic')
        % Using SST of [001], [110], [111]
        [qsort, isort] = sort(abs(qeqv(1, :, i)));
        tmp = qeqv(:, isort, i);
        tmp = tmp(:, reverseCols);
        tmp1 = single(tmp);

        % if sum(tmp1(1, :) < 0)
        %     fprintf('\nMINUS POINT %f\n', tmp1(:, tmp1(1, :) < 0))
        % end

        plusTriangle = find((tmp1(2, :) >= tmp1(3, :)) & (tmp1(3, :) >= tmp1(4, :)) & (tmp1(4, :) >= 0));
        invrTriangle = find((tmp1(2, :) <= tmp1(3, :)) & (tmp1(3, :) <= tmp1(4, :)) & (tmp1(4, :) <= 0));

        % Added this if for Mackenzie cell...
        if isempty(plusTriangle)
            fprintf('\nMINUS POINT %f\n', tmp(:, min(invrTriangle)))
        end
        % newq(:, i) = tmp(:, min(union(plusTriangle, invrTriangle)));
        newq(:, i) = tmp(:, min(plusTriangle));

        %% MORE CHECKS
        if newq(1, i) < cos(deg2rad(62.8)/2)
            disp('hi')
        end
        %%
    end
end

% the Mapped q...
q = newq;
