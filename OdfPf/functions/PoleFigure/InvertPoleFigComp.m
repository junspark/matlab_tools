function varargout = InvertPoleFigComp(odfPfMatrix, pfData, sphL2ipMatrix, FrL2ipMatrix, FrH1sipMatrix, varargin)
% InvertPoleFigComp - This function inverts pole figure data to
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


ztol = 1e-4;    % lower bound zero cutoff

numFrNp = size(FrL2ipMatrix, 1);
frNpWts = sum(FrL2ipMatrix, 1);

poleFigs = pfData.data;

% Optimization options, see "Optimization Toolbox: quadprog"
opts = optimset(...
    'Display', 'off', ...
    'MaxIter', 250);

%%% For use with max entropy test
%%%opts = optimset('Display', 'iter', 'GradObj', 'on', 'DerivativeCheck', 'on');

% ------- defaults -------
emptySphL2ip = 0;

userDisplayFlag = 1;

userDH = 0;
userVMat = 0;

userODFflag = 0;
x0 = ones(numFrNp, 1);

userLBflag = 1;
lb = zeros(numFrNp, 1) + ztol;

userLambda = 0;
lambda = 1;

cell_sphl2ipmatrix = 0;
cell_odfPfMatrix = 0;
cell_pfData = 0;
% ------- defaults -------

% Process options
if nargin > 5
    for i_arg = 1:2:length(varargin)
        if strcmp(varargin{i_arg}, 'InitODF')
            userODFflag = 1;
            x0 = varargin{i_arg + 1};
        elseif  strcmp(varargin{i_arg}, 'nonneg')
            if strcmp(varargin{i_arg + 1}, 'off')
                userLBflag = 0;
                lb = [];
            end
        elseif  strcmp(varargin{i_arg}, 'lb')
            lb = varargin{i_arg + 1};
        elseif  strcmp(varargin{i_arg}, 'Display')
            if strcmp(varargin{i_arg + 1}, 'off')
                userDisplayFlag = 0;
            end
        elseif strcmp(varargin{i_arg}, 'H1Weight')
            userLambda = 0;
            lambda =  varargin{i_arg + 1};
        elseif strcmp(varargin{i_arg}, 'DiscHarm')
            userDH = 1;
            numHarm =  varargin{i_arg + 1};
        elseif strcmp(varargin{i_arg},  'EigMatrix')
            userDH = 1;
            userVMat = 1;
            V =  varargin{i_arg + 1};
            numHarm = size(V, 2);
            if size(V, 1) ~= numFrNp
                error('The discrete harmonic matrix has the wrong number of rows')
            end
        end
    end
end

if userDisplayFlag
    fprintf('\nBeginning pole figure inversion for %s \n**************\n', inputname(2));
end

% Concatenate odfpf matrices into M
if iscell(odfPfMatrix)
    cell_odfPfMatrix = 1;
    numPfs = length(odfPfMatrix);
    if userDisplayFlag
        fprintf('\nConcatenating %d odfpf matrices from %s...', numPfs, inputname(1));
    end
    M = [];
    for i_pf = 1:numPfs
        M = cat(1, M, odfPfMatrix{i_pf});
        numPfNp(i_pf) = size(odfPfMatrix{i_pf}, 1);
    end
    if userDisplayFlag
        fprintf('done\n');
    end
else
    numPfs = 1;
    numPfNp(1) = size(odfPfMatrix, 1);
    M = odfPfMatrix;
end

% Arrange sphere L2IP matrix in block-diagonal L
if iscell(sphL2ipMatrix)
    cell_sphl2ipmatrix = 1;
    if length(sphL2ipMatrix) ~= numPfs
        error('mismatch between sphere L2IP matrix and pole figure data');
    else
        if userDisplayFlag
            fprintf('\nArranging sphere L2IP matricies from %s into block-diagonal form...', inputname(3));
        end
        blSize = zeros(numPfs, 1);
        for i_pf = 1:numPfs
            blSize(i_pf) = size(sphL2ipMatrix{i_pf}, 1);
        end
        L = zeros(sum(blSize));
        ii = 1;
        jj = blSize(1);
        L(ii:jj, ii:jj) = sphL2ipMatrix{1};
        for i_pf = 2:numPfs
            ii = jj + 1;
            jj = jj + blSize(i_pf);
            L(ii:jj, ii:jj) = sphL2ipMatrix{i_pf};
        end
        if userDisplayFlag
            fprintf('done\n');
        end
    end
elseif ~isempty(sphL2ipMatrix) & ~iscell(sphL2ipMatrix) & numPfs > 1
    if userDisplayFlag
        fprintf('\nArranging sphere L2IP matrix from %s into block-diagonal form...', inputname(3));
    end
    L = repmat(full(sphL2ipMatrix), [1, 1, numPfs]);
    L = SparseOfMatArray(L);
    if userDisplayFlag
        fprintf('done\n');
    end
elseif isempty(sphL2ipMatrix)
    if userDisplayFlag
        fprintf('\nRead empty sphere L2IP matrix, handling discrete data\n');
    end
    emptySphL2ip = 1;
else
    L = sphL2ipMatrix;
end

% Concatenate pole figure data into P
if iscell(odfPfMatrix) & iscell(poleFigs)
    cell_pfData = 1;
    if length(poleFigs) ~= numPfs
        error('mismatch between odf -> pf matrices and pole figure data');
    else
        if userDisplayFlag
            fprintf('\nConcatenating Pole Figure data from %s...', inputname(2));
        end
        P = cat(1, poleFigs{:});
        if userDisplayFlag
            fprintf('done\n');
        end
    end
elseif ~iscell(odfPfMatrix) & iscell(poleFigs)
    error('mismatch between odf -> pf matrix and pole figure data');
elseif ~iscell(odfPfMatrix) & ~iscell(poleFigs)
    P = poleFigs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- Computation starts here ------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form objective function for L2 minimization of pole figure data
% (SEE: quadprog.m)
if userDH
    if ~userVMat
        [V, D] = eigs(FrH1sipMatrix, FrL2ipMatrix, numHarm, 'sm');
        [D, Dsort] = sort(diag(D));
        V = V(:, Dsort);
    end
    
    if emptySphL2ip
        H = V'*(M'*M + lambda*FrH1sipMatrix)*V; % for quadprog
        H = 0.5*(H + H');		% force symmetry
        f = -1.0*V'*M'*P;
    else
        H = V'*(M'*L*M + lambda*FrH1sipMatrix)*V; % for quadprog
        H = 0.5*(H + H');		% force symmetry
        f = -1.0*V'*M'*L'*P;
    end
else
    if emptySphL2ip
        H = M'*M + lambda*FrH1sipMatrix; % for quadprog
        f = -1.0*M'*P;			%
    else
        H = M'*L*M + lambda*FrH1sipMatrix;
        f = -1.0*M'*L'*P;
    end
end

% Initial guess and equivalence constraint (for normalization)
if userODFflag
    if userDisplayFlag
        fprintf('\nRead user-input initial ODF estimate\n')
    end
end

% Equality constraint
% **NOTE: for MUD normalization (mean value of unity)
if userDH
    if userLBflag
        A  = -1.0*V;
        b  = -1.0*lb;
    else
        A = [];
        b = [];
    end
    Aeq = frNpWts*V/sum(frNpWts);
    beq = 1;
    x0 = [1/mean(V(:, 1));zeros(numHarm - 1, 1)];
    lb = [];
else
    A  = [];
    b  = [];
    Aeq = frNpWts/sum(frNpWts);
    beq = 1;
end

% Perform optimization
if userDisplayFlag
    fprintf('\nRunning quadratic program for L2 minimization over pole figures\n\n')
    if userLBflag
        fprintf('\t**nonnegativity contraint active\n\n')
    else
        fprintf('\t**nonnegativity contraint inactive\n\n')
    end
end
%tic
[odf, objfun, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, [], x0, opts);
%toc
% EXITFLAG values
%       1  QUADPROG converged with a solution X.
%       3  Change in objective function value smaller than the specified tolerance.
%       4  Local minimizer found.
%       0  Maximum number of iterations exceeded.
%      -2  No feasible point found.
%      -3  Problem is unbounded.
%      -4  Current search direction is not a direction of descent; no further
%           progress can be made.
%      -7  Magnitude of search direction became too small; no further progress can
%           be made.
if exitflag > 0
    if userDisplayFlag
        if exitflag == 1
            mystring = 'QUADPROG converged';
        elseif exitflag == 3
            mystring = 'Change in objective function value smaller than the specified tolerance';
        elseif exitflag == 4
            mystring = 'Local minimizer found.';
        end
        disp(['Optimization terminated successfully: ', mystring]);
    end
elseif exitflag <= 0
    %% debugging %% disp(['we have a problem...', num2str(exitflag)])
    if exitflag == 0
        opts = optimset('Display', 'off', 'MaxIter', 500);
        [odf, objfun, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, [], x0, opts);
        if exitflag == 0
            warning('Maximum number of iterations exceeded in quadprog')
        end
    elseif exitflag == -2
        error('No feasible solution found, exiting')
    elseif exitflag == -3
        error('Problem is unbounded, exiting')
    elseif exitflag == -4
        error('Current search direction is not a direction of descent, exiting')
    elseif exitflag == -7
        error('Magnitude of search direction became too small, exiting')
    end
end

% If using discrete Harmonics, convert to ODF
if userDH
    odf = V*odf;
end

% Output handling
if nargout == 1
    varargout{1} = odf;
elseif nargout == 2
    varargout{1} = odf;
    repfs.hkls = pfData.hkls;
    for i = 1:numPfs
        repfs.data{i} = odfPfMatrix{i}*odf;
    end
    varargout{2} = repfs;
end
