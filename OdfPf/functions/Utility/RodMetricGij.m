function gij = RodMetricGij(refpts)
% RODMETRICGIJ - Compute components of the covariant metric tensor in
% Rodrigues parameter space at a list of reference points.
%   
% USAGE
%     gij = MetricGij(refpts)
%
% INPUTS
%     1) refpts is an 3 x n real array, an array of rodrigues vectors
%
% OUTPUTS
%     1) gij is an 3 x 3 x n real array, the covariant metric tensor
%        components (inner products of the tangent vectors) at each of the
%        n reference points.
%
% NOTES
%     *) The expression used in this function is taken from 
%        Kumar, A., "Modelling crystallographic texture evolution with
%        finite elements over neo-Eulerian orientation spaces", 
%        Comput.  Methods  Appl.  Mech.  Engrg.  153  (1998)  259-302
%
%     *) This routine gives a slightly more accurate metric than one
%        derived from inner products with the tangent vectors output form
%        RodDifferential.  Most likely due to roundoff.
%
% SEE ALSO
%     MetricGij, SphDifferential, RodDifferential
%

npt = size(refpts, 2);
%
ta2    = sqrt(dot(refpts, refpts, 1));
phiby2 = atan(ta2);

ca2 = cos(phiby2);
sa2 = sin(phiby2);

zindex = find(ta2 == 0);

if isempty(zindex)
  axvec = refpts./repmat(ta2, [3, 1]);
else
  ta2(zindex) = 1;
  axvec = refpts./repmat(ta2, [3, 1]);
  axvec(:, zindex) = repmat([1 0 0]', [1, length(zindex)]);
end

ninj = RankOneMatrix(axvec);

T1 = reshape(repmat(reshape(repmat(ca2.^2,       [3, 1]), [1, 3*npt]), [3, 1]), [3, 3, npt]);
T2 = reshape(repmat(reshape(repmat(ca2.^2.*sa2.^2,  [3, 1]), [1, 3*npt]), [3, 1]), [3, 3, npt]);

gij = T1.*repmat(eye(3), [1, 1, npt]) - T2.*ninj;

