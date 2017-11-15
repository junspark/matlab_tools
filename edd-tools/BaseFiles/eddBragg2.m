function  energy = eddBragg2(tth,dhkl)
%
%  calculate Energy of selected hkl (d-spacing)
%
%  p(1) = TTH (degree);
%
% AC@ANL, 4/26/2016

%const = constants;
const.hc = 12.398419057638671;

energy = const.hc/2./dhkl*(1./sind(tth'/2));



