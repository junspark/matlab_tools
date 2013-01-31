function nPrms = nPrmsPkFcn(PkFcn)
if strcmp(PkFcn, 'pkGaussian')
    nPrms   = 3;
elseif strcmp(PkFcn, 'pkLorentzian')
    nPrms   = 3;
elseif strcmp(PkFcn, 'pkpseudoVoigt')
    nPrms   = 4;
elseif strcmp(PkFcn, 'pkpseudoVoigtAdv.m')
    nPrms   = 5;
end