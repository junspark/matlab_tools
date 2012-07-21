function odfpf = SetupOdfPf(pfhkls, mesh, pts, fiberdiv, bsize, varargin)
% SETUPODFPF - make odf -> pf matrices for given set(s) of
% {hkl}, measurement points and odf.
%
% USAGE
%     1) odfpf = SetupOdfPf(pfhkls, mesh, pts, fiberdiv, bsize, options)
%
%
% INPUTS
%     1) pfhkls is 3 x n or 4 x n, the array of Miller or
%         Miller-Bravais indices (plane normals) for the n desired
%         pole figures.
%     2) mesh is a ExtendendMeshStructure over a fundamental region in
%         Rodrigues parameter space containing m independent nodal
%         points.
%     3) pts is 3 x p, the array of pole figure evaluation points
%         (on the unit sphere, usually sphere mesh nodal points)
%         Alternatively, pts may be a 1 x n cell array of p separate
%         coordinate arrays if the coverage varies for each PF
%     4) fiberdiv is 1 x 1, the number of (equispaced) evaluation
%         points for the trapezoid integration of each of the n x p
%         fibers.
%     5) bsize is 1 x 1, the number of pts to evaluate fibers for
%         at once.  This helps to conserve memory for large systems
%         (i.e. n x p > 1e4).
%     6) options is a list of keyword-value pairs:
%         'InversePF',  {'off'}|'on'
%         'RightSym',   'none'|'triclinic'|'orthorhombic'|'tetragonal'|'hexagonal'|{'cubic'}
%         'LeftSym',    {'none'}|'triclinic'|'orthorhombic'
%         'LattPrms',   {ones}|double array
%
% OUTPUTS
%     1) odfpf is a 1 x 3 structure with the following components:
%
%         aspect: the c/a aspect ratio for the unit cell.
%         hkls:   the pfhkls array, for reference.
%         matx:   1 x n cell array of p x m arrays that represent
%                 the centrosymmetric pole figure projection
%                 operators between mesh and pts for each of the n
%                 pole directions.
%
% NOTES
%     Becuase this subroutine operates under the assumption that
%     the orientation space is represented by a mesh structure over
%     the rodrigues parameter space, only the folowing Laue groups
%     are handled:
%
%         orthorhombic:  mmm
%     tetragonal(high):  4/mmm
%       trigonal(high): -3/m
%      hexagonal(high):  6/mmm
%         **cubic(low):  m3
%          cubic(high):  m3m
%
%     **Currently there is no mesh for the cubic(low) symmetry
%     group over rodrigues space.  For the Laue groups
%
%           triclinic: -1
%          monoclinic:  2/m
%     tetragonal(low):  4/m
%       trigonal(low): -3
%      hexagonal(low):  6/m
%
%     the FE/rodrigues parameterization
%     framework cannot be used as the corresponding fundamental
%     regions are unbounded.
%
%     The quaternion parameterization can be used in conjunction
%     with an FE mesh on the projective hypersphere (P3) to handle
%     all Laue groups.  This procedure is under development.
%
% SEE ALSO: ExtendedMeshStructure, FiberOfPoint

% Defaults
userAxiSymFlag = 0;
userSymAxis = [0 0 1]';
cellPts = 0;

lsym = [1 0 0 0]';			% Triclinic default
rsym = CubSymmetries;			% Cubic default
rightSymTag = 'cubic';			%
latticePrms = [];			% initialize lattice parms

%% nrmFac = sum(sum(mesh.l2ip))/2/pi;	% for Odf/Pf normalization*

% Process optional args
if nargin > 5
    for i_arg = 1:2:length(varargin)
        if strcmp(varargin{i_arg}, 'RightSym')
            if strcmp(varargin{i_arg + 1}, 'orthorhombic')
                userRightSymFlag = 1;
                rsym = OrtSymmetries;
                rightSymTag = varargin{i_arg + 1};
            elseif strcmp(varargin{i_arg + 1}, 'trigonal')
                userRightSymFlag = 1;
                rsym = TrigSymmetries;
                rightSymTag = varargin{i_arg + 1};
            elseif strcmp(varargin{i_arg + 1}, 'tetragonal')
                userRightSymFlag = 1;
                rsym = TetSymmetries;
                rightSymTag = varargin{i_arg + 1};
            elseif strcmp(varargin{i_arg + 1}, 'hexagonal')
                userRightSymFlag = 1;
                rsym = HexSymmetries;
                rightSymTag = varargin{i_arg + 1};
            elseif strcmp(varargin{i_arg + 1}, 'cubic')
                userRightSymFlag = 1;
                rsym = CubSymmetries;
                rightSymTag = varargin{i_arg + 1};
            else
                fprintf('\nunrecongnized right symmetry group, using default\n');
                varargin{i_arg}     = CubSymmetries;
                varargin{i_arg + 1} = 'cubic';
            end
        elseif strcmp(varargin{i_arg}, 'LeftSym')
            if strcmp(varargin{i_arg + 1}, 'axisymmetry')
                userAxiSymFlag = 1;
            elseif strcmp(varargin{i_arg + 1}, 'orthorhombic')
                userLeftSymFlag = 1;
                lsym = OrtSymmetries;
                leftSymTag = varargin{i_arg + 1};
            else
                fprintf('\nunrecongnized left symmetry group, using default\n');
            end
        elseif strcmp(varargin{i_arg}, 'symmAxis')
            userSymAxis = varargin{i_arg + 1};
        elseif strcmp(varargin{i_arg}, 'LattPrms')
            latticePrms = varargin{i_arg + 1};
        end
    end
end

npf = size(pfhkls, 2);

% aspect ratio(s) from the lattice parameters
if isempty(latticePrms)
    aspect = 1;
else
    aspect = latticePrms./repmat(latticePrms(1), [1, length(latticePrms)]);
    aspect = aspect(2:end);
end

MBflag = size(pfhkls, 1);
if MBflag == 4
    hkls = MillerBravaisToNormal(pfhkls, aspect);
else
    hkls = pfhkls;
    hkls = MillerIndexToNormal(pfhkls, aspect);
end

fprintf('\nSetupOdfPf: Making ODF-->PF matrices\n')

odfpf.aspect = aspect;

if ~iscell(pts)
    tmpts = {pts};
    clear pts
else
    tmpts = pts;
    cellPts = 1;
    clear pts
end

pind = 1;
for i = 1:npf

    odfpf.hkls(:, i) = pfhkls(:, i);

    hkl = hkls(:, i);

    % first check to see if we need to calculate both +h and -h fibers...
    hsym  = SymHKL(hkl, rsym, 0);
    nhsym = SymHKL(hkl, rsym, 1);

    if size(hsym, 2) == size(nhsym, 2)
        % here +h == -h under the symmetry group
        fprintf('\nProcessing: {%s}, +h == -h', num2str(pfhkls(:, i)'));

        if ~userAxiSymFlag
            odfpf.matx{i} = BuildOdfPfMatrix(hkl, mesh, tmpts{pind}, fiberdiv, bsize, varargin{:});
        else
            odfpf.matx{i} = BuildOdfPfMatrixSph(hkl, mesh, rsym, tmpts{pind}, fiberdiv, bsize, varargin{:});
        end

    elseif size(hsym, 2) ~= size(nhsym, 2)
        % here +h ~= -h under the symmetry group
        fprintf('\nProcessing: {%s}, +h ~= -h', num2str(pfhkls(:, i)'));

        if ~userAxiSymFlag
            plus_odfpf  = BuildOdfPfMatrix(hkl, mesh, tmpts{pind}, fiberdiv, bsize, varargin{:});
            minus_odfpf = BuildOdfPfMatrix(-1.0*hkl, mesh, tmpts{pind}, fiberdiv, bsize, varargin{:});
        else
            plus_odfpf  = BuildOdfPfMatrixSph(hkl, mesh, rsym, tmpts{pind}, fiberdiv, bsize, varargin{:});
            minus_odfpf = BuildOdfPfMatrixSph(-1.0*hkl, mesh, rsym, tmpts{pind}, fiberdiv, bsize, varargin{:});
        end

        % enforce centrosymmetry
        odfpf.matx{i} = 0.5*(plus_odfpf + minus_odfpf);

    end
    if cellPts
        pind = pind + 1;
    end
end

fprintf('\n\nFinished SetupOdfPf, have a nice day!\n')
