function interpMatx = InterpMatrixOfSphMesh(mesh, refpts)
% INTERPMATRIXOFSPHMESH - Evaluate the interpolation matrix for evaluating
%  a function on a sphere mesh a a generic set of coordinate.
%   
% USAGE
%     1) interpMatx = InterpMatrixOfSphMesh(mesh, refpts)
%
% INPUTS
%     1) mesh is a MeshStructure on the unit sphere in dimension d
%         containing l elements, m nodes.
%     2) refpts is the d x n list of n evaluation points.
%
% OUTPUTS
%     1) interpMatx is the (sparse) n x m interpolation matrix
%         that takes the m nodal points of the mesh to the n reference 
%         points.
%
% SEE ALSO EvalSphMeshFunc ProjSphMeshFunc
%
% *) This code could be modified to handle vector fields as well.
%
con = mesh.con;
crd = mesh.crd;

nel = size(con, 2);
dim = size(crd, 1);
nnp = size(crd, 2);
nrf = size(refpts, 2);

%fprintf('\nRead %d points for evaluation\n');

% Make volumetric mesh by adding the origin
newcon = [1 + con;ones(1, nel)];
newcrd = [zeros(dim, 1), crd];

% Move points to tet interiors
refpts = 0.5*UnitVector(refpts);

% Execute search with builtin...
%fprintf('\nPerforming Delaunay search...');
[triList, baryCrds] = tsearchn(newcrd', newcon', refpts');
%fprintf('done\n');
% 
% try to invert if NaN (for hemispheres!)
% *) this could be optimized...to only re-search for the missed pts
if sum(isnan(triList))
    indexer = isnan(triList);
    refpts(:, indexer) = -1.0*refpts(:, indexer);
    [triList, baryCrds] = tsearchn(newcrd', newcon', refpts');
end

baryCrds = baryCrds(:, 1:end-1);
nrm = repmat(sum(baryCrds, 2), [1, dim]);

baryCrds = baryCrds./nrm;
baryCrds = baryCrds';

rowIndex = [1:nrf];
rowIndex = repmat(rowIndex, [dim, 1]);
rowIndex = rowIndex(:);
colIndex = con(:, triList);
colIndex = colIndex(:);

sfun = baryCrds(:);

%fprintf('\nBuild interpolation matrix...');
interpMatx = sparse(rowIndex, colIndex, sfun, nrf, nnp);
%fprintf('done\n');
