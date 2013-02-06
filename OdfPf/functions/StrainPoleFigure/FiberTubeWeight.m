function w = FiberTubeWeight(deta, dome)
% clear all
% close all
% clc
% 
% deta    = 5/2;
% dome    = 10/2;

pts = FiberTube([1 0 0]', deta, dome);
% plot3(pts(1,:), pts(2,:), pts(3,:), 'k.')
% axis equal
% hold on

sph.crd = pts;
sph.con = [1 2 3; 1 3 4; 1 4 5; 1 5 6; 1 6 7; 1 7 8; 1 8 9; 1 9 2]';
sph.eqv = [];
qrule2d = LoadQuadrature('qr_trid06p12');
[~, sph.l2ip]   = SphGQRule(sph, qrule2d);
w   = sum(sph.l2ip)./sum(sum(sph.l2ip));