function [iguess, varargout] = CalculateSPFInversionIG(frmesh, dh, strain, C)
SODF_ISOSTRAIN  = zeros(6,frmesh.numind);

for i = 1:1:frmesh.numind
    R   = RMatOfQuat(QuatOfRod(frmesh.crd(:,i)));
    T   = VectorizedCOBMatrix(R);
    
    C_STAR              = T*C*T';
    SODF_ISOSTRAIN(:,i) = C_STAR*strain;
end

SHIG_11 = GetSHCoeffs(SODF_ISOSTRAIN(1,:), frmesh, dh);
SHIG_22 = GetSHCoeffs(SODF_ISOSTRAIN(2,:), frmesh, dh);
SHIG_33 = GetSHCoeffs(SODF_ISOSTRAIN(3,:), frmesh, dh);
SHIG_12 = GetSHCoeffs(SODF_ISOSTRAIN(4,:), frmesh, dh);
SHIG_13 = GetSHCoeffs(SODF_ISOSTRAIN(5,:), frmesh, dh);
SHIG_23 = GetSHCoeffs(SODF_ISOSTRAIN(6,:), frmesh, dh);

iguess  = [SHIG_11; SHIG_22; SHIG_33; SHIG_12; SHIG_13; SHIG_23]';
iguess  = iguess(:);

nout    = max(nargout,1);
if nout > 1
    thresh1 = 0.5;
    thresh2 = 0.05*max(abs(iguess));
    
    ub  = iguess + abs(iguess)*thresh1;
    lb  = iguess - abs(iguess)*thresh1;
    
    iub = abs(iguess) < thresh2;
    ilb = abs(iguess) < thresh2;
    
    ub(iub) = +thresh2;
    lb(ilb) = -thresh2;
    
    varargout{1}    = lb;
    varargout{2}    = ub;
end
