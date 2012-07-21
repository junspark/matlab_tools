%
%  ODF From Pole Figures.
%
%  This script illustrates the methods available for
%  constructing an ODF from pole figure data.
%
%
%-------------------- User Input
%
wsname = 'wspfi';      % workspace name
h1tol  = 1.0e-1;
%
%-------------------- Execution
%
%  Load workspace for fundamental region.
%
addpath('../build-workspaces');
%
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname);
%
% 1) Find best pole figure.
%
% X=QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0)
%
sphl2 = ws.sphmesh.l2ip;
%
H = ws.pfmats(1).odfpf'*sphl2*ws.pfmats(1).odfpf + ...
    ws.pfmats(2).odfpf'*sphl2*ws.pfmats(2).odfpf;
%
f = -1*(...
    ws.pfmats(1).odfpf'*sphl2*ws.pfs(1).data + ...
    ws.pfmats(2).odfpf'*sphl2*ws.pfs(2).data  ...
    );
%
A   = [];
b   = [];
Aeq = full(sum(ws.frmesh.l2ip));
beq = sum(Aeq);
LB  = zeros(1, ws.frmesh.numind);
UB  = [];
X0  = ones(ws.frmesh.numind, 1);
%
disp('running first optimization');
odf1 = quadprog(H,f,A,b,Aeq,beq,LB,UB,X0);
bestpf1 = ws.pfmats(1).odfpf * odf1;
bestpf2 = ws.pfmats(2).odfpf * odf1;
%
% 1) Find best ODF.
%
% X=QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0)
%
H = ws.frmesh.h1form;
f = zeros(ws.frmesh.numind, 1);
%
A1 = [...
    ws.pfmats(1).odfpf; ...
    ws.pfmats(2).odfpf  ...
    ];
b1  = A1*odf1;
tol = h1tol*ones(size(b1));
%
A  = [A1; -A1];
b  = [b1 + tol; -b1 + tol];
%
%  Aeq and beq as before.
%
X0 = odf1;
%
disp('running second optimization');
odf = quadprog(H,f,A,b,Aeq,beq,LB,UB,X0);
finalpf1 = ws.pfmats(1).odfpf * odf;
finalpf2 = ws.pfmats(2).odfpf * odf;
%
%  Create DX output files.
%
Ndata = {'odf1', odf1, 'odf', odf};
ExportDX('pfs-odf', ws.frmesh, Ndata);
%
Ndata = {'raw-110', ws.pfs(1).data, 'raw-100', ws.pfs(2).data, ...
	'best-110', bestpf1, 'best-100', bestpf2, ...
	'final-110', finalpf1, 'final-100', finalpf2};
%
ExportDX('pfs-odf-pfs', ws.sphmesh, Ndata);
save
