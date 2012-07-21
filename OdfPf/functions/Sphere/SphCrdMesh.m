function smesh = SphCrdMesh(ntheta, nphi, varargin)
% SphCrdMesh - Generate a hemisphere mesh based on spherical coordinates.
%   
%   USAGE:
%
%   smesh = SphCrdMesh(ntheta, nphi)
%
%   INPUT:
%
%   ntheta is a positive integer,
%          the number of subdivisions in theta
%   nphi   is a positive integer,
%          the number of subdivisions in phi
%
%   These arguments can be followed by a list of
%   parameter/value pairs which specify keyword options.
%   Available options are listed below with default values
%   shown in brackets.
%
%   'MakeL2ip'     {'on'}|'off'
%                  computes the L2 inner product matrix and 
%                  adds it to the mesh structure as a field .l2ip
%   'QRule'        string  {'qr_trid06p12'}
%                  name of a quadrature rule for triangles, to be used
%                  in building the l2ip matrix
%
%   OUTPUT:
%
%   smesh is a MeshStructure,
%         on the hemisphere (H^2)
%
%   NOTES:
%
%   * No equivalence array is produced.
%

%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'MakeL2ip',  'on', ...
    'QRule',     'qr_trid06p12'
       };
%
opts = OptArgs(optcell, varargin);
%
%-------------------- Execution
%
%
tdiv = (0:ntheta)*2*pi/ntheta;
pdiv = (0:nphi)*pi/(2*nphi);
%
[phi, theta] = meshgrid(pdiv, tdiv);
npts = (ntheta+1)*(nphi+1);
%
thetaphi = [reshape(theta, 1, npts) ; reshape(phi  , 1, npts)];
%
xyz = XYZOfThetaPhi(thetaphi);
%
%  Now clean up coordinates to remove duplicates
%  due to spherical coordinate parameterization.
%
nt1 = ntheta + 1;
np1 = nphi + 1;
%
leftedge  = 1:nt1:1+nt1*nphi;
rightedge = leftedge + ntheta;
%
SeeNode = 1:npts;
SeeNode(rightedge) = leftedge;
SeeNode(1:nt1) = 1;
%
UseThese = (SeeNode >= 1:npts);
nreduced = sum(UseThese);
%
scrd = xyz(:, UseThese);    % coordinates
%
NewNode = (1:npts);
Masters = NewNode(UseThese);
NewNode(Masters) = 1:nreduced;
%
%  Now set up connectivity and eliminate degenerate elements.
%
top = (1+nt1*nphi):npts;
OldNodes = 1:npts;
OldNodes(rightedge) = 0;
OldNodes(top) = 0;
NodeOne = OldNodes(OldNodes > 0);
%
tcon1 = [NodeOne; NodeOne + 1; NodeOne + nt1 + 1];
tcon2 = [NodeOne; NodeOne + nt1 + 1; NodeOne + nt1];
tmpcon = NewNode(SeeNode([tcon1 tcon2]));
%
Eq12 = (tmpcon(2,:) - tmpcon(1,:) == 0);
Eq13 = (tmpcon(3,:) - tmpcon(1,:) == 0);
Eq23 = (tmpcon(2,:) - tmpcon(3,:) == 0);
%
Degenerate = max([Eq12 ; Eq13; Eq23]);
%
scon = tmpcon(:, ~Degenerate);
%
smesh = MeshStructure(scrd, scon, []);
%
%  Compute L2ip if desired.
%
if OnOrOff(opts.MakeL2ip)
  [tmp, smesh.l2ip] = SphGQRule(smesh, LoadQuadrature(opts.QRule));
end
