function pseudo_strain = funcEDDInstrObjective(p)
% funcEDDInstrObjective - EDD instrument function used to determine the
% instrument parameters
%
%   USAGE:
%
%   pseudo_strain = funcEDDInstrObjective(p)
%
%   INPUT:
%
%   p
%       instrument parameters; p(1) is 2theta in degrees and p(2:end) are
%       the polynominal coefficients to convert channels to photon energy.
% 
%   OUTPUT:
%
%   pseudo_strain
%       strain computed from the known d-spacing to the measured d-spacing.
%

data    = load('obj_data.mat');
data    = data.data;

TOA                     = p(1);
ChToEnergyConversion    = p(2:end);

x           = data.x;
y           = data.y;
peaks2use   = data.peaks2use;
d_hkl       = data.d_hkl;

lambda_hkl0 = 2.*d_hkl*sind(TOA/2);
E_hkl0      = Angstrom2keV(lambda_hkl0);
E_grid      = polyval(ChToEnergyConversion, x);

for i = 1:1:length(peaks2use)
    E0  = E_hkl0(peaks2use(i));

    idx1    = find(E_grid < (E0 - 1.5));
    idx2    = find(E_grid < (E0 + 1.5));
    idx1    = idx1(end);
    idx2    = idx2(end);
    
    xdata   = E_grid(idx1:idx2);
    ydata   = y(idx1:idx2);
    
    p0  = [ ...
        max(ydata); ...
        2.5; ...
        0.5; ...
        E0; ...
        0; ...
        0; ...
        ];
    
    pLB = [0; 0; 0; E0-3; -inf; -inf];
    pUB = [inf; inf; 1; E0+3; inf; inf];
    
    p   = lsqcurvefit(@pfunc, p0, xdata, ydata, pLB, pUB);
    
    A(i)    = p(1);
    
    yfit0   = pfunc(p0, xdata);
    yfit    = pfunc(p, xdata);
        
    E_fit(i)    = p(4);

end
pseudo_strain   = E_fit./E_hkl0(peaks2use) - 1;