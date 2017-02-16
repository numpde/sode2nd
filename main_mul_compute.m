function main_mul_compute
% function main_mul_compute
%
%   Convergence study for 2nd moment / covariance
%   with multiplicative noise
%
%   Writes the results to main_mul_results.mat
%
%   Reference: https://arxiv.org/abs/1611.02164
%    
%   Requires the ppfem package:
%   https://bitbucket.org/numpde/ppfem/downloads
%   Extract the zip to folder DIR and use addpath('DIR/lib/');

%   Author: R. Andreev, Aug 2016

    % Helpers
    iif = @(C, varargin) varargin{2-logical(C)};
    Vec = @(X) X(:);
    
    % ppfem path
    ifexists = @(p) iif(exist(p, 'dir'), p, '');
    addpath(ifexists('./../../lib'));
    
    
    %%%%%%%%%%%%%%%%%%%%
    % Problem setup
    %%%%%%%%%%%%%%%%%%%%

    % SODE:
    % dX + lambda X = rho X dB  with  X(0) = X0

    % Final time
    problemdata.T = 2;
    % Initial datum
    problemdata.X0 = 1;
    % Spectral parameter
    problemdata.lambda = 3;
    % Volatility
    problemdata.rho = 0.5 * sqrt(2*problemdata.lambda);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Convergence
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Method encoding:
    % ----------------
    % +(p)      = CN* family of degree p for the test space F 
    % -(p)      = iE* family of degree p for F with solution postprocessing
    % -(p+1000)	= iE*/Q = like "-(p)" but with Delta Q
    % -(p+2000)	= iE*/box = like "-(p)" but with "quadrature" for Delta
    
    figure;
    
    styles = {'bo-', 'b*-', 'rs-', 'rx-', 'gd-', 'mo-'};
    labels = {'CN_2^*', 'CN_2^*(2)', 'iE_2^*', 'iE_2^*(2)', 'iE_2^*/Q', 'iE_2^*/box'};
    methods = [+(1), +(2), -(1), -(2), -(1 + 1000), -(1 + 2000)];

    % Debugging options:
%     styles = {'mo-'}; labels = {'CN*'}; methods = [(1)];
%     styles = {'mo-'}; labels = {'CN*(2)'}; methods = [(2)];
%     styles = {'mo-'}; labels = {'CN*(3)'}; methods = [(3)];
%     styles = {'mo-'}; labels = {'P iE_2^*'}; methods = [-(1)];
%     styles = {'mo-'}; labels = {'P iE_2^*(2)'}; methods = [-(2)];
%     styles = {'mo-'}; labels = {'P iE_2^*/Q'}; methods = [-(1 + 1000)];
%     styles = {'mo-'}; labels = {'P iE_2^*/\box'}; methods = [-(1 + 2000)];

    % Output filename
    filename = 'main_mul_results';

    % Postprocessing option
    do_postprocessing = true;
    if (do_postprocessing)
        methods = methods + (sign(methods) * 1e4);
    end
    disp(['Postprocessing: ' num2str(do_postprocessing)]);

    % Number of temporal mesh elements
    NN = 2 .^ (3:9);
    for method = methods
        KK = NaN(size(NN));
        EE = NaN(size(NN));
        
        disp(['Computing: ' labels{method == methods}]);
        for N = NN
            mesht = linspace(0, problemdata.T, N+1);
            KK(N == NN) = max(diff(mesht));
            
            [RES, REF] = compute(problemdata, mesht, method);
            
            % 2nd moment
            f = @(t) Vec(REF.M2(t,t))';
            g = @(t) Vec(ppfem_eval(RES.diag_PPE, t) * RES.diag_PP_M2)';
            
%             % Covariance
%             f = @(t) Vec(REF.CovX(t,t))';
%             g = @(t) Vec(ppfem_eval(RES.diag_QE, t) * RES.diag_Q_CovX)';
            
            % Compute trace error in L1
            e = integral(( @(t) abs(f(t) - g(t)) ), 0, problemdata.T, 'RelTol', 1e-1);
            
            disp(['Error: ' num2str(e)]);
            EE(N == NN) = e;
    
            % Show the diagonal
            tt = linspace(0, problemdata.T, 1e3);
            semilogy(tt, f(tt), '--', tt, g(tt), '-');
            title('Diagonal'); xlabel('t'); pause(0.1);
            
            clear mesht RES REF f g e tt
        end
            
        KE{method == methods}.KK = KK;
        KE{method == methods}.EE = EE;
        
        disp('--');
        
        save([filename '.mat'], '-mat');
    end
end


function [RES, REF] = compute(problemdata, mesht, method)
% function [RES, REF] = compute(problemdata, mesht, method)

    %%%%%%%%%%%%%%%%%%%%
    % Helpers
    %%%%%%%%%%%%%%%%%%%%
    
    Vec = @(X) X(:);
    Row = @(X) X(:)';
    Mat = @(X) reshape(X, sqrt(numel(X)) * [1, 1]);
    Out = @(x, y) Vec(x) * Row(y); % "outer" product

    %%%%%%%%%%%%%%%%%%%%
    % Problem parameters
    %%%%%%%%%%%%%%%%%%%%

    % Deterministic initial value 
    X0 = problemdata.X0; 
    % Stochastic ODE parameters
    lambda = problemdata.lambda;
    rho = problemdata.rho;

    % Exact expected value
    REF.EX   = @(t) X0 * exp(-lambda*t);
    % Exact 2nd moment
    REF.M2   = @(s,t) (X0^2 * exp(-lambda*(s+t)) .* (exp(rho^2 * min(s,t))));
    % Exact covariance
    REF.CovX = @(s,t) (X0^2 * exp(-lambda*(s+t)) .* (exp(rho^2 * min(s,t)) - 1));

    %%%%%%%%%%%%%%%%%%%%
    % FEM parameters
    %%%%%%%%%%%%%%%%%%%%

    % % Refine towards the origin
    % mesht = sort(unique([mesht, mesht(2) * (2.^(-(1:3)))]));

    % Piecewise polynomial degree P >= 1 of the test space
    P = mod(abs(method), 1000);

    % Construct trial space
    if (method == 0)
        error('Unknown method = 0');
    elseif (method > 0)
        % CN* family
        E = ppfem_construct(mesht, P-1, 'L');
    elseif (method < 0)
        % iE* family and variants
        E = ppfem_construct(mesht, P, 'R*');
    
        % % Scale such that <e, f> evaluates f (optional)
        % E = ppfem_recombine(ppfem_eval(E, E.x(1:end-1), '+') / ppfem_gramian(E), E);
    end

    % Construct test space: IP^P space with v(T) = 0
    F = ppfem_construct(mesht, P, 'P');
    F = ppfem_homconstraint(F, ppfem_eval(F, mesht(end)));
    
    % Mass matrices
    withmass = @(V) setfield(V, 'M', ppfem_gramian(V));
    E = withmass(E);
    F = withmass(F);
    dF = withmass(ppfem_transform(F, @fnder));

    % Additional temporal FEM matrices
    F_x_E = ppfem_gramian(F, E);
    dF_x_E = ppfem_gramian(dF, E);

    % Evaluation operator at t = 0 (row vector)
    f0 = ppfem_eval(F, [0]);
    
    % Evolution operator (acting on the unknown, not the test function)
    b = -dF_x_E  +  lambda * F_x_E;
    
    do_postprocessing = (abs(method) >= 1e4);
    method = sign(method) * mod(abs(method), 1e4);
    
    % Preprocessed 1d solution space
    QE = ppfem_construct(mesht, P-1, 'L');
    % Preprocessing operator = 
    % projection onto piecewise polynomials
    Q = ppfem_gramian(QE) \ ppfem_gramian(QE, E);
    

    %%%%%%%%%%%%%%%%%%%%
    % COMPUTE E[X]
    %%%%%%%%%%%%%%%%%%%%

    % Compute first moment EX = IE[X]
    % The variational formulation is
    % Find EX in X_lambda s.t.
    %     int <EX, (-d_t + lambda) v> dt = <IE[X0], v(0)>
    % for all  v in Y_lambda

    EX = b \ (X0 * Vec(f0));

    %%%%%%%%%%%%%%%%%%%%
    % COMPUTE Cov[X]
    %%%%%%%%%%%%%%%%%%%%

    % Compute covariance CovX = Cov[X]
    % The variational formulation is
    % Find CovX in X_lambda^(2) s.t.
    %     int int <CovX, (-d_s + lambda) (-d_t + lambda) v(s,t)> ds dt - int <rho^(2) CovX(t,t), v(t,t)> dt 
    %     =
    %     <Cov[X0], v(0,0)> + int <(rho E[X])^(2)(t,t), v(t,t)> dt
    % for all v in Y_lambda^(2)

    % Norm on the trial space
    ME = lambda * E.M;
    % Norm on the test space
    MF = (1/lambda) * dF.M  +  (lambda) * F.M  +  Out(f0, f0);
    
    [diag_E, E_tt] = vtt(E);
    [diag_F, F_tt] = vtt(F);
    [diag_QE, QE_tt] = vtt(QE);
    
    %%%%%%%%%%%%%%%%%%%%
    % Operators
    %%%%%%%%%%%%%%%%%%%%
    
    % Case-by-case Delta (transpose) operator
    if ((0 < abs(method)) && (abs(method) < 1000))
        % CN* / iE* families exact Delta (transpose)
        Int0T = ppfem_gramian(diag_E, diag_F);
        Delta = @(C) Mat(F_tt' * Int0T' * E_tt * Vec(C));
        Delat = @(Y) Mat(E_tt' * Int0T  * F_tt * Vec(Y));
        clear Int0T;
        
    elseif ((1000 < -method) && (-method < 2000))
        % iE*/Q family: Include projection Q inside Delta
        Int0T_Q = ppfem_gramian(diag_QE, diag_F);
        Delta = @(C) Mat(F_tt' * Int0T_Q' * QE_tt * Vec(Q * C * Q'));
        Delat = @(Y) Q' * Mat(QE_tt' * Int0T_Q  * F_tt * Vec(Y)) * Q;
        clear Int0T_Q;
        
    elseif ((2000 < -method) && (-method < 3000))
        % iE*/box
        assert(P == 1, 'iE*/box undefined unless P = 1');
        
        fn = (ppfem_gramian(E) / ppfem_eval(E, E.x(1:end-1), '+'))' * ppfem_eval(F, mesht(1:end-1));
        Delta = @(C) fn' * diag(diag(C) ./ diff(mesht)') * fn;
        Delat = @(Y) diag(diag(fn * Y * fn') ./ diff(mesht)');
        clear fn;
    
    else
        error('Unknown method');
    end
    
    % Operator "B - rho^2 Delta" and its transpose
    B  = @(C) (b * C * b'  -  rho^2 * Delta(C));
    Bt = @(Y) (b' * Y * b  -  rho^2 * Delat(Y));
    
    % Preconditioners
    iM2 = @(C) (ME \ C / ME');
    iN2 = @(Y) (MF \ Y / MF');
    
    % Solver
    Solve = @(ell) Mat(pcg(@(C)Vec(Bt(iN2(B(Mat(C))))), Vec(Bt(iN2(ell))), 1e-10, 1000, @(C)Vec(iM2(Mat(C)))));
    
    %%%%%%%%%%%%%%%%%%%%
    % Solve
    %%%%%%%%%%%%%%%%%%%%
    
    %%% 2nd moment %%%
    
    % Right hand side for the 2nd moment
    ell = X0^2 * Out(f0, f0);
    
    % Solve B(CovX) = b with preconditioning
    M2 = Solve(ell);
    
    %%% Covariance %%%
    
    % Right hand side for the covariance
    ell = rho^2 * Delta(Out(EX, EX));
    % (Not clear if this is the correct ell in all cases)
    
    % Solve B(CovX) = b with preconditioning
    CovX = Solve(ell);
    
    
    %%%%%%%%%%%%%%%%%%%%
    % Collect output
    %%%%%%%%%%%%%%%%%%%%
    
    if do_postprocessing
        PPE = QE;
        PP = Q;
        diag_PPE = diag_QE;
        PPE_tt = QE_tt;
    else
        PPE = E;
        PP = 1;
        diag_PPE = diag_E;
        PPE_tt = E_tt;
    end
    
    RES.problemdata = problemdata;
    RES.mesht = mesht;
    RES.method = method;
    RES.do_postprocessing = do_postprocessing;
    %
    RES.E = E;
    RES.EX = EX;
    RES.M2 = M2;
    RES.CovX = CovX;
    %
    RES.PPE = PPE;
    RES.PP_EX = PP * EX;
    RES.PP_M2 = PP * M2 * PP';
    RES.PP_CovX = PP * CovX * PP';
    %
    RES.diag_PPE = diag_PPE;
    RES.diag_PP_M2 = PPE_tt * Vec(RES.PP_M2);
    RES.diag_PP_CovX = PPE_tt * Vec(RES.PP_CovX);
end


function [diag_V, Vtt] = vtt(V)
% function [diag_V, Vtt] = vtt(V)
%
%   Vtt * C(:) returns the coefficients of the diagonal of C
%   wrt the space diag_V

    % Helpers
    Row = @(X) X(:)';
    getfields = @(s, name) [s(:).(name)];
    
    % Procedure: 
    %   1. V --> locally orthogonal Legendre polynomial space U
    %   2. U x U --> diag(U x U) --> diag_V

    maxP = (max(getfields([V.pp{:}], 'order')) - 1);
    U = ppfem_construct(V.x, maxP, 'L');
    
    max2P = 2 * maxP;
    diag_V = ppfem_construct(U.x, max2P, 'L');
    
    % The mapping (2.) on the reference element [0,1]
    LLL = reshape(intI_LiLjLk([0,1], maxP), [(maxP+1)^2, max2P+1])';
    
    % Assemble the mapping (2.) into Utt
    Utt = sparse(diag_V.dim, U.dim^2);
    % Interval sizes
    k = diff(U.x);
    for n = 1 : length(k)
        % Indices of basis functions on n-th interval
        bn = @(P) (1 + (0:P) + (n-1) * (P+1));
        % ii = "output" indices of basis functions of diag_V space
        ii = bn(max2P);
        % jj = "input" indices of operator
        jj = Row(sub2ind([1,1]*U.dim, meshgrid(bn(maxP)), meshgrid(bn(maxP))'));
        % Transfer from the reference element
        Utt(ii, jj) = LLL / sqrt(k(n));
    end
    
    U.M = 1; % Mass matrix of U = Id
    V2U = U.M \ ppfem_gramian(U, V);
    Vtt = Utt * kron(V2U, V2U);
end

function LLL = intI_LiLjLk(I, P)
%function LLL = intI_LiLjLk(I, P)
%
%   Integral of triple products of 
%   orthonormal Legendre polynomials 
%   on the interval I = [a, b]
%
%   The degrees are
%       0 <= i <= P
%       0 <= j <= P
%       0 <= k <= 2 P
%
%   The size of the resulting tensor LLL is
%       (P+1) x (P+1) x ((2*P)+1)
%
%   Thus:
%       Li Lj = sum_k LLL(i,j,k) Lk
%
%   Might be useful:
%       reshape(LLL, [(P+1)^2, (2*P)+1])
%   and
%       sqrt(diff(I)) * intI_LiLjLk(I, P) == intI_LiLjLk([0,1], P)

%   R. Andreev, 2016-08-18

    pp = getfield(ppfem_construct(I, 2*P, 'L'), 'pp');
    LLL = zeros(P+1, P+1, (2*P)+1);
    for i = (0 : P) + 1
        for j = (0 : P) + 1
            for k = (0 : 2*P) + 1
                LLL(i, j, k) = pp0integral(pp0cmb(pp{i}, '*', pp0cmb(pp{j}, '*', pp{k})));
            end
        end
    end
    clean = @(A) (A .* (abs(A) > 1e-10 * max(abs(A(:)))));
    LLL = clean(LLL);
end
