%
%  Evaluate spherical harmonics on meshes.
%
%  * Evaluate both on S^2 (dim=3) and S^3 (dim=4)
%  * Degree(S^3) = 2* degree(S^2)
%
%-------------------- * Input
%
DEG = 4;  % degree on S^2
WSName = 'c4xhs15x';
DXFileName = sprintf('harmonics-deg-%d-%s.dx', DEG, WSName);
%
SaveWS   = 0;
SaveName = 'Workspace';
%
%-------------------- * Execution
%
load(WSName);
%
[cf3, pb3] = PolyHarmonicBasis(3, DEG);
[cf4, pb4] = PolyHarmonicBasis(4, 2*DEG);
%
%  Now, evaluate.
%
%  S^2
%
fprintf('Working on S^2\n')
npts  = size(ws.sphmesh.crd, 2);
nhrm  = size(cf3, 2);
%
harmonics.s2 = zeros(npts, nhrm);
for i=1:nhrm
  fprintf('evaluating %d of %d\n', i, nhrm);
  harmonics.s2(:, i) = PolyBasisEval(pb3, cf3(:, i), ws.sphmesh.crd);
end
%
%  S^3
%
fprintf('Working on S^3\n')
npts  = ws.frmesh.numind;
nhrm  = size(cf4, 2);
%
qpts   = QuatOfRod(ws.frmesh.crd(:, 1:ws.frmesh.numind));
harmonics.s3 = zeros(npts, nhrm);
for i=1:nhrm
  fprintf('evaluating %d of %d\n', i, nhrm);
  harmonics.s3(:, i) = SymmHarmonic(pb4, cf4(:, i), qpts, ws.frmesh.symmetries);
end
%
%  Find basis.
%
if (max(abs(harmonics.s3(:))) < 1.0e-8)
  harmonics.s3sym = zeros(npts, 1);
else
  harmonics.s3sym = orth(harmonics.s3);
end
%
%WriteDXHarmonics(DXFileName, ws, harmonics);
%
%  Save workspace.
%
if SaveWS
   save(SaveName); 
end
%-------------------- * Clean up
%
%clear ...

