function WriteDXHarmonics(fname, ws, h)
% WriteDXHarmonics - Write IsaiahMesh in DX format.
%   
%  fname is a string
%           the name of the output file (include .dx suffix)
%  ws is a Workspace structure
%  h is a structure (with .s2 and .s3 giving harmonics)
%
%  * 
%

%
%dxf = DXFile('Name', fname);
%
%-------------------- Fundamental Region Field
%
NData = cell(2, size(h.s3sym, 2));
for i=1:size(h.s3sym, 2)
  NData{1, i} = sprintf('FR Harmonic %d', i);
  NData{2, i} = h.s3sym(:, i)';
end
%NData = reshape(NData, [1, 2*size(NData, 2)]);
frfld = ExportDXField([fname, '-fr'], ws.frmesh, NData);
%
%-------------------- Sphere Field
%
NData = cell(2, size(h.s2, 2));
for i=1:size(h.s2, 2)
  NData{1, i} = sprintf('S2 Harmonic %d', i);
  NData{2, i} = h.s2(:, i)';
end
%NData = reshape(NData, [1, 2*size(NData, 2)]);
sphfld = ExportDXField([fname, '-sph'], ws.sphmesh, NData);
return
%
