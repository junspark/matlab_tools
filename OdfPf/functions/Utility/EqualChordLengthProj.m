function pcrd = EqualChordLengthProj(crd, basis)
% EQUALCHORDLENGTHPROJ - Generalization of the 2-d Equal-area
% projection for meshes on the projective n-sphere
% 
%  USAGE:
%          pcrd = EqualChordLengthProj(crd, basis)


dim  = size(crd, 1);
npts = size(crd, 2);

refbasis = eye(dim);
pdir = refbasis(:, dim);

if (nargin == 1)
    pcrd = crd;   
else
    if size(basis, 1) ~= dim
        error('inconsistent basis')
    end
    pcrd = basis'*crd;     % components in basis
end

chords = pcrd - repmat(pdir, [1, npts]);
cnorms = NormVecArray(chords);
cnorms(cnorms == 0) = 1;

pchords = chords(1:dim - 1, :);
pcnorms = NormVecArray(pchords);
pcnorms(pcnorms == 0) = 1;

pcrd = repmat(cnorms, [dim - 1, 1]).*chords(1:dim - 1, :)./repmat(pcnorms, [dim - 1, 1]);
