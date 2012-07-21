function outMesh = FixConnectivity(inMesh)
% FIXCONNECTIVITY - Find and repair and inside-out elements (det(J) < 0) in
% a connectivity array
%
% USAGE
%     outMesh = FixConnectivity(inMesh)
%
% INPUTS
%     1) inMesh is a MeshStructure with a simplex element connectivity
%
% OUTPUTS
%     1) outMesh is a MeshStructure identical to inMesh, with any
%        inside-out elements repaired.  
%
% NOTES
%     *) An inside-out element is defined as one with a negative Jacobian
%        determinant.  MATLAB's builtin delaunay tesselation algorithm can
%        produce inside-out simplices.  It is recommended to run this
%        algorithm on any mesh obtained by delaunayn.m.
%

outMesh = inMesh;

Jac = Jacobian(inMesh);

negJac = find(Jac < 0);

if isempty(negJac)
  fprintf('\nFixConnectivity found no inside-out elements.\n');
else
  len = length(negJac);
  
  fprintf('\nFixConnectivity found %d inside-out elements\n\trepairing...', ...
    len); 
  for i = 1:len
    outMesh.con(:, negJac(i)) = [inMesh.con(1:end-2, negJac(i));
                                 inMesh.con(end, negJac(i));
                                 inMesh.con(end-1, negJac(i));
                                ];
  end
  fprintf('done\n');
end