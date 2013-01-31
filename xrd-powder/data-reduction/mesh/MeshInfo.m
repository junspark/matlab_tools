function MeshInfo(mesh, sym, tol)
% MeshInfo - Verify and print properties of mesh.
%   
%   USAGE:
%
%   outinfo = MeshInfo(mesh)
%   outinfo = MeshInfo(mesh, sym)
%   outinfo = MeshInfo(mesh, sym, tol)
%
%   INPUT:
%
%   mesh is a MeshStructure
%   sym  is 4 x n,  (optional)
%        the symmetry group in quaternions; this argument can 
%        be omitted if the equivalence array is empty or if you 
%        do not desire to verify it
%   tol  is a scalar, (optional, default = 1.0e-7)
%        the tolerance used in verifying and generating 
%        equivalences 
%
%   OUTPUT:  none
%
%   This function prints various information to the screen.
%
crd = mesh.crd;
con = mesh.con;
eqv = mesh.eqv;
%
[mcrd, ncrd] = size(crd);
[mcon, ncon] = size(con);
[meqv, neqv] = size(eqv);
%
if (nargin < 3)
  tol = 1.0e-7;
end
%
indent = '    ';
fprintf(1, '--- Nodal Points\n');
fprintf(1, '%s%37s%d\n', ...
	indent, 'dimension of nodal points:  ', mcrd, ...
	indent, 'number of nodal points:  '   , ncrd, ...
	indent, 'number of independent nodal points:  ', (ncrd-neqv) ...
	);
%
fprintf(1, '\n');
fprintf(1, '--- Equivalences\n');
%
if (isempty(eqv))
  %
  fprintf(1, '%s%27s\n', indent, ...
	  'Equivalence array is empty.');
  %
else
  %
  if (meqv ~= 2) 
    error('First dimension of equivalence array is not 2.')
  end
  %
  if (nargin < 2)
    %
    fprintf(1, '%s%s%s\n%s%s\n', ...
	  indent, 'Equivalence array is nonempty, ', ...
	          'but no symmetry group is given.', ...
	  indent, 'Cannot verify equivalences.'  ...
	    );
    nsym = 0;
    %
  else
    %
    if (isempty(sym))
      %
      fprintf(1, '%s%s\n%s%s\n', ...
	      indent, 'Empty symmetry group given.', ...
	      indent, 'Not verifying equivalences.'  ...
	      );
      nsym = 0;
      %
    else
      nsym = size(sym, 2);
    end
    %
  end
  %
  if ((neqv > 0)&(nsym > 0) )
    %
    fprintf(1, '%s%s%d\n', indent, ...
	    'number of equivalences: ', neqv);
    
    fprintf(1, '\n');
    fprintf(1, '%s... Verifying equivalences ... ', indent);
    noneqv = VerifyEqv(mesh, sym, tol);
    fprintf(1, 'done\n');
    fprintf(1, '%s    (tolerance used:  %e)\n', indent, tol);
    %
    if (isempty(noneqv))
      fprintf(1, '%s--- All equivalences are valid.\n', indent);
    else
      fprintf(1, '%s%s\n%s%s%d\n%s%s\n', ...
	      indent, '*** Invalid equivalences exist.', ...
	      indent, '*** number of invalid equivalences:  ', (length(noneqv)), ...
	      indent, '*** see "outinfo.noneqv" for the list of bad columns');
      outinfo.noneqv = noneqv;  % output 
    end
    
    fprintf(1, '\n');
    fprintf(1, '%s... Generating equivalences independently ... ', indent);
    nmesh = ReduceMesh(mesh, sym, tol);
    fprintf(1, 'done\n');
    fprintf(1, '%s    (tolerance used:  %e)\n', indent, tol);

    neqvnew = size(nmesh.eqv, 2);
    if (neqvnew == neqv)
      fprintf(1, '%s--- Number of equivalences agrees with input mesh..\n', ...
	      indent);
    else
      fprintf(1, '%s%s\n%s%s%d\n%s%s%d\n', ...
	      indent, '*** Check equivalence arrray.', ...
	      indent, '***    number of computed equivalences:  ', neqvnew, ...
	      indent, '***    number of original equivalences:  ', neqv ...
	      );
    end
  end
end
%
fprintf(1, '\n');
fprintf(1, '--- Connectivity\n');
%
fprintf(1, '%s%33s%d\n', ...
	indent, 'dimension of reference element:  ', (mcon -  1), ...
	indent, 'number of elements:  ', ncon ...
	);
%
fprintf(1, '\n');
fprintf(1, '%sTesting for mismatched interior faces.\n', indent);
fprintf(1, '%s... Generating list of faces ... ', indent);
[faces, multiplicity] = MeshFaces(con);
fprintf(1, 'done\n');
nmult1 = sum(multiplicity == 1);
nmult2 = sum(multiplicity == 2);
nmult3 = sum(multiplicity >= 3);
if (nmult3 > 0)
  fprintf(1, '%s%s\n%s%s\n', ...
	indent, '*** Check your connectivity.',...
	indent, '*** There are element faces with multiplicity greater than two.' ...
	);
else
    fprintf(1, '%s%40s%d\n%s%40s%d\n', ...
	indent, '    number of matched element faces:  ',       nmult2, ...
	indent, '    number of unmatched element faces:  ', nmult1  ...
	);
end
%
scon   = faces(:, (multiplicity == 1));
[faces, multiplicity] = MeshFaces(scon);
nmult1 = sum(multiplicity == 1);
if (nmult1 == 0)
  fprintf(1, '%s--- Found no unmatched interior faces.\n', indent);
else
  fprintf(1, '%s*** Check your mesh.\n', indent);
  fprintf(1, '%s*** There are unmatched interior faces.\n', indent);
end
%
if(neqv > 0)
fprintf(1, '\n');
fprintf(1, '%sTesting symmetries of surface elements.\n', indent);
fprintf(1, '%s... Applying equivalences to connectivity ...', ...
	indent);
smesh = MeshStructure(crd, scon, eqv);
sconr = ReduceConnectivity(smesh);
fprintf(1, 'done\n');
%
nsel  = size(sconr, 2);
scon1 = [ones(1, nsel); 1+sconr];
[f1, m1] = MeshFaces(scon1);
select = (min(f1, [], 1) >  1);
m1 = m1(select);
%
nmult1 = sum(m1 == 1);
nmult2 = sum(m1 == 2);
nmult3 = sum(m1 >= 3);
%
if (nmult3 > 0)
  fprintf(1, '%s*** Check your mesh.\n', indent);
  fprintf(1, '%s%s%s\n',...
	  indent, ...
	  '*** There are element faces with ', ...
	  'multiplicity greater than two.' ...
	  );
else
  if (nmult1 > 0)
    fprintf(1, '%s%s\n%s%s%d\n',...
	    indent, '*** There are unmatched faces.', ...
	    indent, '***    number of unmatched faces:  ', nmult1 ...
	    );
  else
    fprintf(1, '%s%s\n',...
	    indent, '--- There are no unmatched faces.' ...
	    );
  end
  %
end
end
%
%
