function rmat = RMatOfGRAINDEXU(U)
  % RMATOFGRAINDEXU - A utility routine to reformate the row-wise list of
  % orientation matrix components output from GRAINDEX into a 3 x 3 x n array of
  % rotation matrices.
  %
  % USAGE:
  %
  %      rmat = RMatOfGRAINDEXU(U)
  %
  % INPUTS:
  %
  %      1) U is n x 9, a row-concatenated list of n orientation matrix
  %      components.
  %
  % OUTPUTS:
  %
  %      1) rmat is 3 x 3 x n, an array of n rotation matrices
  %
  % SEE ALSO
  %
  
  numUs = size(U, 1);
  
  rmat = permute(reshape(U', [3 3 numUs]), [2 1 3]);
end
