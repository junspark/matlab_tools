%
%  Make pole figures.
%
%-------------------- User Input
%
bground = 0.5;  % background (% of uniform)
centers = [
    0.25    -0.5
    0.25     0.3
    0.25     0.0
    ];
%
stdevs  = [0.3, 0.2];     % standard deviations
%
%-------------------- Execution
%
%  Load a workspace.
%
addpath('../build-workspaces'); % use \ for windows
load wshex
%
%  Create sample ODF.
%
allpts = wshex.frmesh.crd;
nindep = wshex.frmesh.numind;
pts    = allpts(:, 1:nindep);  % independent nodes 
%
symm  = wshex.frmesh.symmetries;
odf   = zeros(1, nindep);
cenfr = ToFundamentalRegion(QuatOfRod(centers), symm);
%
for i=1:size(centers, 2)
  rg  = RodGaussian(cenfr(:, i), pts, stdevs(i), symm);
  odf = odf + rg ./ MeanValue(rg, wshex.frmesh.l2ip);
end
%
odf  = odf ./ MeanValue(odf, wshex.frmesh.l2ip);
%
unif = ones(1, nindep);
unif = unif ./ MeanValue(unif, wshex.frmesh.l2ip);
%
odf = bground*unif + (1-bground)*odf;
%
%  Plot the odf.
%
fropts = {'Symmetries', 'hexagonal', 'ShowMesh', 'on'};
PlotFR(wshex.frmesh, odf, fropts{:});
%
%  Make and display the pole figure.
%
pf100 = wshex.pfmats(1).odfpf*odf(:);
%
plotopts = {'UpperHemisphere', 'on'};
PlotSphere(wshex.sphmesh, pf100, plotopts{:});
%
pf001 = wshex.pfmats(2).odfpf*odf(:);
%
PlotSphere(wshex.sphmesh, pf001, plotopts{:});
%
pfdata = {pf100, eye(3), pf001, eye(3)};
%
PlotPF2d(wshex.sphmesh, pfdata);
%
%  Create DX output files.
%
Ndata = {'odf', odf};
ExportDX('hexodf', wshex.frmesh, Ndata);
%
Ndata = {'pf100', pf100, 'pf001', pf001};
ExportDX('hexpfs', wshex.sphmesh, Ndata);
%
save
