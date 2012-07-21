function cvec = ToAxisymmetricFR(quat, qsym, symmAxis)

nsym = size(qsym, 2);
ncrd = size(quat, 2);

Rmat = RMatOfQuat(quat);
Rsym = SparseOfMatArray(RMatOfQuat(qsym));

ztol  = 1e-8;

for i = 1:ncrd
  % crd(:, i) = Rmat(:, :, i)'*symmAxis;
  % for j = 1:nsym
  %   csym{i}(:, j) = Rsym(:, :, j)*crd(:, i);
  % end
  %% if strcmp(leftSymTag, 'cubic')
  %% Using SST of [100], [110], [111]
  % tmpcrd  = csym{i};
  tmpcrd = Rmat(:, :, i)'*symmAxis;
  tmpcrd = reshape(Rsym*repmat(tmpcrd, [nsym, 1]), [3, nsym]);

  tmpcrd(find(abs(tmpcrd) < ztol)) = 0;

  tmpcrd1 = single(tmpcrd);

  SST1 = find((tmpcrd1(1, :) >= tmpcrd1(2, :)) & (tmpcrd1(2, :) >= tmpcrd1(3, :)) & (tmpcrd1(3, :) >= 0));
  SST2 = find((tmpcrd1(2, :) >= tmpcrd1(1, :)) & (tmpcrd1(1, :) >= tmpcrd1(3, :)) & (tmpcrd1(3, :) >= 0));
  %
  if isempty(SST1) & isempty(SST2)
    disp('wait!  We have a roundoff problem')
  end
  cvec(:, i) = tmpcrd(:, min(union(SST1, SST2)));
  %
end