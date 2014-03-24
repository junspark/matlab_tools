function mesh = MeshStructureXRD(crd, con, numel, qrule)
% MESHSTRUCTUREXRD - Create mesh structure from mesh data for XRD image
%   
%   mesh = MeshStructure(crd, con, numel, eqv)
%
%   crd is e x n, the nodal point array 
%   con is d x m, the connectivity
%   numel is scalar, the number of elements
%   qrule is 3 x k the equivalence array
%
%   mesh is a mesh-structure, which consists of the fields
%        above; 
%
mesh = struct('crd', crd, ...
    'con', con, ...
    'numel', numel, ...
    'qrule', qrule);
