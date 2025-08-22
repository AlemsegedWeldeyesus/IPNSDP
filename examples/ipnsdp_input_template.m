function ipnsdp_input_template()
%IPNSDP_INPUT_TEMPLATE  Example template for defining an NSDP in ipnsdp.
%
% This script sets up a generic Nonlinear Semidefinite Program (NSDP) 
% of the form:
%
%   minimize    f(X, u)
%   subject to
%       sum_{k=1}^{m_X} <A11_i^{(k)}, X_k> + A12_i * u = b1_i,   i = 1,...,n1
%       sum_{k=1}^{m_X} <A21_j^{(k)}, X_k> + A22_j * u <= b2_j,  j = 1,...,n2
%       c1(X, u) = 0
%       c2(X, u) <= 0
%       A_s(X, u) >= 0   (PSD),   s = 1,...,n_psd
%       u_lower <= u <= u_upper
%       rho_k_lower * I_{d_k} <= X_k <= rho_k_upper * I_{d_k},   k = 1,...,m_X
%
% where 
%   - X_k in S^{d_k},  k = 1,...,m_X, are symmetric matrix variables
%   - u in R^n is a vector of scalar variables
%   - <A, X> = trace(A' * X) denotes the Frobenius inner product
%
%% Problem Definition
% All fields are listed with empty values. Comments explain intended usage.

prob.name         = '';         % Name of the problem (optional)
prob.problem_data = struct();   % Struct to store user data (optional)

% Matrix variables X
prob.nX      = [];              % m_X: number of matrix variables X_k. 
                                % Set to 0 if none.
prob.dimX    = [];              % [d_1,...,d_m_X]: sizes d_k of X_k. 
                                % Set to [] if none.

% Bounds on Matrix variables X
prob.lbX    = [];               % [rho_1_lower,...,rho_m_X_lower]
                                % Set to [] if all unbounded. 
                                % Use -Inf for selectively unbounded entries.
prob.ubX    = [];               % [rho_1_upper,...,rho_m_X_upper]
                                % Set to [] if all unbounded.     
                                % Use +Inf for selectively unbounded entries.

% Scalar variables u
prob.nu      = [];              % n: dimension of vector u. 
                                % Set to 0 if none.

% Bounds on scalar variables u
prob.lbu     = [];              % [u_lower_1;...;u_lower_n]
                                % Set to [] if all unbounded.
                                % Use -Inf for selectively unbounded entries.
prob.ubu     = [];              % [u_upper_1;...;u_upper_n]
                                % Set to [] if all unbounded.
                                % Use +Inf for selectively unbounded entries.

% nPSD constraints
prob.nPSDcon = [];              % n_psd: number of PSD constraints A_s(X,u) >= 0.
                                % Set to 0 if none.

% Linear constraints
% Note: You may embed these into nonlinear constraints if easier.
prob.A11 = [];       % [vec(A11_1_1)^T, ..., vec(A11_1_mX)^T;
                     %  ...;
                     %  vec(A11_n1_1)^T, ..., vec(A11_n1_mX)^T]
                     % Coefficient matrix for [vec(X_1); ...; vec(X_mX)]
                     % in linear equality constraints. [] if none.
prob.A12 = [];       % Coefficient matrix for scalar variables u in 
                     % linear equalities. [] if none.
prob.b1  = [];       % RHS vector of linear equality constraints. 
                     % [] if none.
prob.A21 = [];       % [vec(A21_1_1)^T, ..., vec(A21_1_mX)^T;
                     %  ...;
                     %  vec(A21_n2_1)^T, ..., vec(A21_n2_mX)^T]
                     % Coefficient matrix for [vec(X_1); ...; vec(X_mX)]
                     % in linear inequality constraints. [] if none.
prob.A22 = [];       % Coefficient matrix for scalar variables u in 
                     % linear inequalities. [] if none.
prob.b2  = [];       % RHS vector of linear inequality constraints. 
                     % [] if none.

% Objective function
prob.f_obj = @(x) objective(x, prob);  
                                % User-defined objective function 
% Nonlinear and PSD constraints
prob.c = @(x) nonlinear_constraint(x, prob);  
                     % User-defined nonlinear constraints
                     % Must return [eq, ineq, psd]:
                     % eq   = nonlinear equalities (vector)
                     % ineq = nonlinear inequalities (vector)
                     % psd  = cell array of PSD matrices

% Initial point (optional)
prob.x0.X = {};                 % Cell array of initial guesses for X_k. 
                                % Example: prob.x0.X{1} = eye(d_1), ...,
                                %          prob.x0.X{m_X} = eye(d_m_X)
                                % Defaults are used if not specified.
prob.x0.u = [];                 % Column vector initial guess for u. 
                                % Default is zero if not specified.

%% Algorithmic Options
% Default solver parameters (modifiable by the user)

prob.options.scale            = 'on';   % 'on'/'off': Enable/Disable scaling
prob.options.opt_tol          = 1e-6;   % Optimality tolerance
prob.options.feas_tol         = 1e-7;   % Feasibility tolerance
prob.options.max_iter         = 100;    % Maximum number of iterations
prob.options.min_mu           = 1e-9;   % Minimum barrier parameter
prob.options.sigma            = 0.4;    % Centrality parameter
prob.options.tau              = 0.9;    % Backtracking step size
prob.options.gamma            = 0.99;   % Step length scaling
prob.options.stepMode         = 'joint';% 'joint'/'separate': Step strategy 
prob.options.chord_decomp     = 'yes';  % 'yes'/'no': Chordal decomposition
prob.options.clique_merge     = 'yes';  % 'yes'/'no': Clique merging
prob.options.density          = 0.7;    % Density threshold for decomposition
prob.options.min_block_size   = 30;     % Minimum block size for decomposition
prob.options.rel_block_size   = 0.6;    % Relative block size for decomposition
prob.options.linconCheck      = 'yes';  % 'yes'/'no': LinConst dependency test
prob.options.nonlinconCheckN  = 'yes';  % 'yes'/'no': NonlinConst dependency test (Numeric)
prob.options.nonlinconCheckS  = 'no';   % 'yes'/'no': NonlinConst dependency test (Symbolic)
prob.options.penlpMode        = 'off';  % 'off'/'on': Penalization
prob.options.penlp            = 0;      % Penalty parameter (0 = off)
prob.options.penldeg          = 2;      % 1 or 2: Penalty degree 

%% Solve the problem
prob.solve = @() ipnsdp_solve(prob);
[x_sol, info] = prob.solve();

%% Results (saved to ipnsdp_sol.mat)
x_sol.X              % Cell array of solution X_k matrices or slacks 
x_sol.u              % Optimized scalar variables vector u
x_sol.Z              % Dual vars for PSD constraints X_k and A_s(X,u)
x_sol.lambda1        % Dual variables for linear equalities
x_sol.lambda2        % Dual variables for linear inequalities
x_sol.theta1         % Dual variables for nonlinear equalities
x_sol.theta2         % Dual variables for nonlinear inequalities
x_sol.slinear        % Slack variables for linear inequalities
x_sol.snonlinear     % Slack variables for nonlinear inequalities

info.obj_evals                     % objective function evaluations
info.grad_evals                    % gradient evaluations
info.nonlincon_evals               % nonlinear constraint evaluations
info.nonlincon_jacobi_evals        % nonlinear constraint Jacobian evals
info.Lag_hessian_evals             % Lagrangian Hessian evaluations
info.CPU                           % CPU time used (seconds)
info.iter                          % Number of iterations performed
info.obj_value                     % Final objective function value
info.opt_tol                       % Final optimality tolerance
info.primal_feasibility            % Final primal feasibility measure
info.dual_feasibility              % Final dual feasibility measure
info.complementarity               % Final complementarity gap
info.obj_val_hist                  % History of objective values
info.kkt_optimality                % History of KKT conditions
info.exit                          % Solver exit message

end

%% Objective Function
function f_val = objective(x, prob)
% Computes the objective function value.
%
% Inputs:
%   x   : struct with fields:
%           x.X : cell array of matrix variables
%           x.u : vector of scalar variables
%   prob: problem-specific data struct
%
% Output:
%   f_val: scalar objective value
%
% ... (user-defined implementation)
end

%% Nonlinear and PSD Constraints
function [eq, ineq, psd] = nonlinear_constraint(x, prob)
% Computes nonlinear equality, inequality, and PSD constraints.
%
% Inputs:
%   x.X : cell array {X_1, ..., X_mX}, each X_k in S^{d_k}
%   x.u : vector of scalar variables u in R^n
%   prob: problem data struct
%
% Outputs:
%   eq   : nonlinear equality constraints vector (if none set to [])
%   ineq : nonlinear inequality constraints vector (if none set to [])
%   psd  : psd{1} = A_1(X, u),..., psd{n_psd} = A_n_psd(X, u) 
%          cell array of PSD constraint matrices (if none set to {})
%
% ... (user-defined implementation)
end
