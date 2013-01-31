function [f, Jac] = pkpseudoVoigtAdv(p, x)

p   = p(:);
x   = x(:);

A       = p(1);
Gamma   = p(2);
etaA    = p(3);
etaB    = p(4);
xPeak   = p(5);

eta     = etaA + etaB*x;

pG  = [A Gamma xPeak];
pL  = [A Gamma xPeak];

[fG, JacG]  = pkGaussian(pG, x);
[fL, JacL]  = pkLorentzian(pL, x);

f   = eta.*fG + (1-eta).*fL;

% compute jacobian
dfdA        = eta.*JacG(:,1) + (1-eta).*JacL(:,1);
dfdGamma    = eta.*JacG(:,2) + (1-eta).*JacL(:,2);
dfdetaA     = fG - fL;
dfdetaB     = x.*fG - x.*fL;
dfdxPeak    = eta.*JacG(:,3) + (1-eta).*JacL(:,3);

% arrange in Jacobian matrix
Jac = [...
    dfdA, ...
    dfdGamma, ...
    dfdetaA, ...
    dfdetaB, ...
    dfdxPeak, ...
    ];