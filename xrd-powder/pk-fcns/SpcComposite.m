function [yCalc, Jac] = SpcComposite(p, x, nPks, M, PkType, bckOrder)
nx      = length(x);

PkFcn   = PkType{1};
PkPrms  = PkType{2};

Jac_bck = zeros(nx, bckOrder+1);
Jac_pks = zeros(nx, nPks*PkPrms);

p_bck   = p(1:bckOrder+1);
p_pks   = p(bckOrder+2:end);

% compute y for background
yCalc_bck   = polyval(p_bck, x);
for i = 1:1:bckOrder
    Jac_bck(:,i)    = x.^(bckOrder-i+1);
end
Jac_bck(:,bckOrder+1) = 1;

% Lorentz polarization factor
L	= LPfactor(x./2);

yCalc_pks   = x*0;
for i = 1:1:nPks
    
    i_ini   = PkPrms*(i-1) + 1;
    i_fin   = PkPrms*i;
    p_pk    = p_pks(i_ini:i_fin);
    
    [yCalc_pk, Jac_pk]      = PkFcn(p_pk, x);
    Jac_pks(:,i_ini:i_fin)	= Jac_pk;
        
    yCalc_pk    = yCalc_pk*M(i);
    Jac_pks(:,i_ini:i_fin)  = Jac_pks(:,i_ini:i_fin)*M(i);
    
    yCalc_pks   = yCalc_pks + yCalc_pk;
end
yCalc_pks   = yCalc_pks.*L;

L       = L*ones(1,nPks*PkPrms);
Jac_pks = Jac_pks.*L;

yCalc   = yCalc_pks + yCalc_bck;
Jac     = [Jac_bck Jac_pks];