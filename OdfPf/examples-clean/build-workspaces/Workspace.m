function ws = Workspace(frmesh, sphmesh, varargin)
% Workspace - Make Odf/Pf workspace.
%
%   USAGE:
%
%   ws = Workspace(frmesh, sphmesh, 'param', 'value')
%
%   INPUT:
%
%   frmesh  is a MeshStructure,
%           on a fundamental region; it is expected to
%           have an additional structure component, 'symmetries'
%   sphmesh is a MeshStructure,
%           on the sphere
%
%   OUTPUT:
%
%   ws is a structure
%
ws.frmesh  = frmesh;
ws.frmesh.numind = size(frmesh.crd, 2) - size(frmesh.eqv, 2);
%
ws.sphmesh = sphmesh;
%
%--------------------Defaults and Options-------------------------------
%
optcell = {...
    'PointsPerFiber',   100,  ...
    'MakePoleFigures', {{}},  ...
    'MakeFRL2IP',     'off',  ...
    'MakeFRH1IP',     'off',  ...
    'MakeSphL2IP',    'off',  ...
    'MakeSphH1IP',    'off',  ...
    'ProgressReport',  'on',  ...
    'QRuleTetrahedra', 'qr_tetd06p24',  ...
    'QRuleTriangles',  'qr_trid06p12',  ...
    'NoOption',          []
    };
%
opts = OptArgs(optcell, varargin);
%
report = OnOrOff(opts.ProgressReport);
%
%-------------------- Execution
%
%  Make requested pole figures.
%
pf_hkls   = opts.MakePoleFigures;
if (~isempty(pf_hkls))
    ws.pfmats = struct('hkl', pf_hkls);
end
%
for i=1:length(pf_hkls)
    %
    hkl   = pf_hkls{i};
    ppf   = opts.PointsPerFiber;
    block = max(floor(5.0e4/ppf), 1); % 50,000 points at once
    %
    if report
        disp(['working on OdfPf matrix:  ', num2str(i)]);
        disp(['   hkl:  ', num2str(hkl(:))']);
    end
    %
    ws.pfmats(i).odfpf = ...
        BuildOdfPfMatrix(hkl, ...
        frmesh, frmesh.symmetries, ...
        sphmesh.crd, ppf, 0, block);
    %
end
%
%  Make various quadratic forms.
%
%  *** Fundamental Region
%
makefrl2 = OnOrOff(opts.MakeFRL2IP);
makefrh1 = OnOrOff(opts.MakeFRH1IP);
%
if (makefrl2 | makefrh1)
    qrule3d = LoadQuadrature(opts.QRuleTetrahedra);
    gqrule  = QRuleGlobal(frmesh, qrule3d, @RodMetric);
    npqp    = NpQpMatrix(frmesh, qrule3d);
    ws.frmesh.qrule = qrule3d;
end
%
if (makefrl2)
    if report
        disp(['working on FR L2IP matrix:  ']);
    end
    ws.frmesh.l2ip = L2IPMatrix(gqrule, npqp);
end
%
if (makefrh1)
    if report
        disp(['working on FR H1IP matrix:  ']);
    end
    ws.frmesh.h1form = H1SIPMatrix(frmesh, qrule3d, ...
        RodDifferential(frmesh, qrule3d.pts));
end
%
%  *** Sphere
%
makesphl2 = OnOrOff(opts.MakeSphL2IP);
makesphh1 = OnOrOff(opts.MakeSphH1IP);
%
if (makesphl2 | makesphh1)
    qrule2d = LoadQuadrature(opts.QRuleTriangles);
    gqrule  = QRuleGlobal(sphmesh, qrule2d, @RodMetric);
    npqp    = NpQpMatrix(sphmesh, qrule2d);
    ws.sphmesh.qrule = qrule2d;
end
%
if (makesphl2)
    if report
        disp(['working on Sphere L2IP matrix:  ']);
    end
    [gqr, ws.sphmesh.l2ip] = SphGQRule(sphmesh, qrule2d);
end
%
if (makesphh1)
    if report
        disp(['working on Sphere H1SIP matrix:  ']);
    end
    ws.sphmesh.h1form = SphH1SIP(sphmesh, qrule2d);
end
%
return
