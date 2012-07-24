function s = C2S(c, varargin)
%%% INPUT C (SHEAR IN TERMS OF SIGXY = C44*EPSXY)
%%% OUTPUT S (SHEAR IN TERMS OF EPSXY = C44*SIGXY)

% default options
optcell = {...
    'Symmetry', 'cubic', ...
    };

% update option
opts	= OptArgs(optcell, varargin);

if strcmpi(opts.Symmetry, 'cubic')
    c11 = c(1);
    c12 = c(2);
    c44 = c(3);     %%% SHEAR MOD C(3) IS IN TERMS OF SIGXY = C(3)*EPSXY
                    %%% NO NEED TO CONVERT c44
    m   = c11*c11 + c11*c12 - 2*c12*c12;
    
    s11 = (c11 + c12)/m;
    s12 = -c12/m;
    s44 = 1/c44;    %%% ALREADY IN EPSXY = s44*SIGXY
    
    s   = [s11; s12; s44];
elseif strcmpi(opts.Symmetry, 'hexagonal')
    c11	= c(1);
    c33 = c(2);
    c12 = c(3);
    c13 = c(4);
    
    c44 = c(5);     %%% SHEAR MOD C(5) IS IN TERMS OF SIGXY = C(3)*EPSXY
    c44 = c44/2;    %%% CONVERT C44 TO SIGXY = C44*GAMXY
    
    m1  = c11 - c12;
    m2  = c33*(c11 + c12) - 2*c13*c13;
    
    s11 = 0.5*(1/m1 + c33/m2);
    s12 = 0.5*(c33/m2 - 1/m1);
    s13 = -c13/m2;
    s33 = (c11 + c12)/m2;
    
    s44 = 1/c44;
    s44 = s44/2;    %%% CONVERT s44 TO EPSXY = s44*SIGXY
    
    s   = [s11; s33; s12; s13; s44];
else
    error('Symmetry not supported')
end