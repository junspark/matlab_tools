function gdotv = GradDotVel_Form(m, qrule)
% GradDotVel_Form.m -- form which integrates gradient against a vector field
%
%   USAGE:
%
%   gdotv = GradDotVel_Form(m, qrule)
%
%   INPUT:
%
%   m     is a MeshStructure
%            the mesh
%   qrule is a QRuleStructure
%            a qudrature rule (local)
%
%   OUTPUT:
%
%   gdotv is n x 3n, sparse
%       it is the matrix which integrates (grad f . V)
%
%   NOTES:
%
%   * Rodrigues metric is used here; this needs to be modified
%     if needed for usual Euclidean geometries
%
numEl  = size(m.con, 2);
numQp  = length(qrule.wts);
numSf  = 4 ; % linear tets
numAll = numEl * numQp * numSf;
%
qruleGl = QRuleGlobal(m, qrule, @RodMetric);
grad    = QuatGradient(m, qrule);
%
%  Build reorientation velocity at nodes/dof, 
%  i.e.  4 (cmps) x  3 dof/node x  4 nodes
%
fprintf('building reorientation velocity field ... ')
q = SpreadRefPts(m, qrule.pts);
%
% now, sf @ qp, want identity matrix at each el/qp/sf
%
sfqp = qrule.pts(:)';
allVecs = repmat(eye(3), [1, 1, numSf, numQp, numEl]) .* ...
	  reshape(repmat(sfqp(:)', [9 1 numEl]), [3, 3, numSf, numQp, numEl]);
quats = reshape(QuatOfRod(qruleGl.pts), [4 numQp numEl]);
fprintf('done\n')
%
%  Put them together
%
fprintf('building integration matrix ... ')
eMats = zeros(numSf, 3*numSf, numEl);
eWts  = reshape(qruleGl.wts, [numQp numEl]);
numDof = 3*numSf;
dof1 = 1:3:numDof; dof2 = 2:3:numDof; dof3 = 3:3:numDof;
clear qruleGl
for i=1:numEl
  % tangent vectors on quaternion space
  eReor  = reshape(allVecs(:, :, :, :, i), [3 numDof numQp]);
  eQuats = repmat(quats(:, :, i),          [1 numDof]);
  % no different % eQuats = permute(repmat(quats(:, :, i), [1 1 numDof]), [1 3 2]);
  eQvecs = QuatReorVel(reshape(eQuats, [4 numDof numQp]), eReor);
  % gradients:  copy to give equal number
  qpGrad = grad(:, :, :, i);
  % take dot product and sum over qp
  edotP = zeros(numSf, 3*numSf);
  wts   = eWts(:, i);
  for j=1:numSf
    for k=1:numDof
      dotP = dot(eQvecs(:, k, :), qpGrad(:, j, :));
      eDotP(j, k) = sum(squeeze(dotP).*wts);
    end
  end
  % insert in elemental matrix list
  eMats(:, :, i) = eDotP;
end
% finally, make a sparse matrix out of it
gi = repmat(m.con, [numDof 1]);
%
dofCon = [3*m.con(:)' - 2; 3*m.con(:)' - 1; 3*m.con(:)'];
gj = repmat(dofCon(:)', [numSf 1]);
%
gMat = sparse(gi, gj, eMats(:));
% now reduce by equivalence
gMat = EqvReduce(gMat', m.eqv);
%
dofEqv = [3*m.eqv-2, 3*m.eqv-1, 3*m.eqv];
gdotv = EqvReduce(gMat', dofEqv);

fprintf('done\n')
