function varargout = InvertPoleFigCompScl(odfPfMatrix, pfData, sphL2ipMatrix, FrL2ipMatrix, FrH1sipMatrix, varargin)
% INVERTPOLEFIGMINL2 - This function inverts pole figure data towhich
% obtain an ODF whose pole figures are the optimal fit of the input
% data (in the L2-norm sense).
%
% USAGE
%     1) [odf, repfs] = InvertPoleFigComp(odfPfMatrix, pfData, sphL2ipMatrix, FrL2ipMatrix, FrH1sipMatrix, options)
%
% INPUTS
%     1) odfPfMatrix is
%     2) pfData
%     3) sphL2ipMatrix
%     4) FrL2ipMatrix
%     5) options are string-value pairs:
%         Initial guess:              'InitODF', [ double ]
%         Non-negativity constraint:  'nonneg', [ 'off' | {'on'} ]
%
% OUTPUTS
%     1) odf
%     2) repfs
%
%   SEE ALSO:  SetupOdfPf, L2IPMatrix, H1SIPMatrix
%

warning off all

fprintf('\nBeginning pole figure inversion for %s \n**************\n', inputname(2));

clear global ODF
global ODF

ztol = 1e-3;    % lower bound zero cutoff

numFrNp = size(FrL2ipMatrix, 1);
frNpWts = sum(FrL2ipMatrix, 1);
numPfs  = length(pfData.data);

% Optimization options, see "Optimization Toolbox: quadprog"
opts = optimset(...
  'Display', 'iter', ...
  'TolFun', ztol, ...
  'TolX', ztol, ... 
  'DiffMinChange', 1.0e-3, ...
  'DiffMaxChange', 1.0e-1, ...
  'MaxIter', 20);
  
  % 'MaxIter', 20, ...

% ------- defaults -------
userInitflag = 0;
x0 = ones(numPfs, 1);

userODF = 0;
if isempty(ODF)
    ODF = ones(numFrNp, 1);
end

userDH = 0;
userVMat = 0;

userLBflag = 1;
%% lb = zeros(numPfs, 1) + ztol;
lb = ztol*x0;
ub = []; % 100*x0;

userLambda = 0;
lambda = 1;
% ------- defaults -------

% Process options
if nargin > 5
    for i_arg = 1:2:length(varargin)
        if strcmp(varargin{i_arg}, 'InitGuess')
            userInitflag = 1;
            x0 = varargin{i_arg + 1};
        elseif strcmp(varargin{i_arg}, 'InitODF')
            userODF = 1;
            ODF =  varargin{i_arg + 1};
        elseif  strcmp(varargin{i_arg}, 'nonneg')
            if strcmp(varargin{i_arg + 1}, 'off')
                userLBflag = 0;
                lb = [];
            end
        elseif strcmp(varargin{i_arg}, 'H1Weight')
            userLambda = 0;
            lambda =  varargin{i_arg + 1};
        elseif strcmp(varargin{i_arg}, 'DiscHarm')
            userDH = 1;
            numHarm =  varargin{i_arg + 1};
        elseif strcmp(varargin{i_arg},  'EigMatrix')
            userVMat = 1;
            userDH = 1;
            V =  varargin{i_arg + 1};
            numHarm = size(V, 2);
        end
    end
end

if userDH
    if ~userVMat
        [V, D] = eigs(FrH1sipMatrix, FrL2ipMatrix, numHarm, 'sm');
        [D, sortD] = sort(diag(D));
        V = V(:, sortD);
    end
    tic
    alfa = fmincon(...
        @PFI_nrmFac_func, x0, ...
        [], [], ...
        [], [], ...
        lb, ub, ...
        [], opts, ...
        odfPfMatrix, pfData, ...
        FrL2ipMatrix, FrH1sipMatrix, ...
        lambda, numHarm, V);
    toc

else
    tic

    alfa = fmincon(...
        @PFI_nrmFac_func, x0, ...
        [], [], ...
        [], [], ...
        lb, ub, ...
        [], opts, ...
        odfPfMatrix, pfData, ...
        FrL2ipMatrix, FrH1sipMatrix, ...
        lambda);
    toc

end

%Assign Solution
odf.fnc = ODF;
odf.scl = alfa;

% Output handling
if nargout == 1
    varargout{1} = odf;
elseif nargout == 2
    repfs.hkls = pfData.hkls;
    for i = 1:numPfs
        repfs.data{i} = odfPfMatrix{i}*ODF;
    end
    varargout{1} = odf;
    varargout{2} = repfs;
end

%% OLD
%     for i = 1:numPfs
%         sclPfData.hkls(:, i) = pfData.hkls(:, i);
%         sclPfData.data{i} = alfa(i)*pfData.data{i};
%     end
%
%     [odfSol, repfs] = InvertPoleFigComp(...
%         odfPfMatrix, sclPfData, [], ...
%         FrL2ipMatrix, FrH1sipMatrix, ...
%         'InitODF', ODF, ...
%         'H1Weight', lambda, ...
%         'DiscHarm', numHarm, ...
%         'EigMatrix', V, ...
%         'Display', 'off');
