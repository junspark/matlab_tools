function [p0, p_LB, p_UB] = GenerateSpcIGuess(x, y, theta_hkl, PkType, bckOrder)
% GenerateSpcIGuess - generates initial guess for a spectrum.
%
%   USAGE:
%
%   p0 = GenerateSpcIGuess(x, y, theta_hkl, PkType, bckOrder)
%
%   INPUT:
%
%   x
%       x coordinates typically in 2 theta, mm, or pixels
%
%   y
%       intensity value corresponding to x
%
%   theta_hkl
%       theta of the peaks present in x
%
%       theta_hkl
%
%   PkType
%       cell array specifying number of peak parameters and peak function
%       type
%
%   bckOrder
%       polynomial order for background
%
%   OUTPUT:
%
%   p0
%       initial guess of the spectrum
%
%   p_LB
%       lower bound for spectrum fit parameters
%
%   p_UB
%       upper bound for spectrum fit parameters

nPks    = length(theta_hkl);

PkPrms  = PkType{2};
PkFcnStr    = char(PkType{1});

% background initial guess
p0_bck  = zeros(bckOrder+1,1);
p0_bck(end) = y(1);
p_bck_LB    = -Inf*ones(bckOrder+1,1);
p_bck_UB    = Inf*ones(bckOrder+1,1);

p0_pks  = zeros(nPks, PkPrms);
p_pks_LB    = -Inf*ones(nPks, PkPrms);
p_pks_UB    = Inf*ones(nPks, PkPrms);

for i = 1:1:nPks
    if strcmp(PkFcnStr, 'pkGaussian')
        A_pk        = .25;
        gamma_pk    = 0.1;
        tth_pk      = theta_hkl(i)*2;

        p0_pks(i,:) = [A_pk gamma_pk tth_pk];
    elseif strcmp(PkFcnStr, 'pkLorentzian')
        A_pk        = 0.5;
        gamma_pk    = 0.03;
        tth_pk      = theta_hkl(i)*2;

        p0_pks(i,:) = [A_pk gamma_pk tth_pk];
    elseif strcmp(PkFcnStr, 'pkpseudoVoigt')
        A_pk        = 0.25;
        gamma_pk    = 0.03;
        eta         = 0.5;
        tth_pk      = theta_hkl(i)*2;

        p0_pks(i,:)     = [A_pk gamma_pk eta tth_pk];
        p_pks_LB(i,:)   = [-Inf -Inf 0 -Inf];
        p_pks_UB(i,:)   = [Inf Inf 1 Inf];
    elseif strcmp(PkFcnStr, 'pkpseudoVoigtAdv')
        A_pk        = 0.25;
        gamma_pk    = 0.03;
        etaA        = 0.5;
        etaB        = 0;
        tth_pk      = theta_hkl(i)*2;

        p0_pks(i,:)     = [A_pk gamma_pk etaA etaB tth_pk];
        p_pks_LB(i,:)   = [-Inf -Inf -Inf -Inf -Inf];
        p_pks_UB(i,:)   = [Inf Inf Inf Inf Inf];
    end
end

p0_pks      = p0_pks';
p_pks_LB    = p_pks_LB';
p_pks_UB    = p_pks_UB';

p0      = [p0_bck(:); p0_pks(:)];
p_LB    = [p_bck_LB(:); p_pks_LB(:)];
p_UB    = [p_bck_UB(:); p_pks_UB(:)];