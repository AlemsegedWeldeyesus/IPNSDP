function mc_pd_sdp()
% MC_PD_SDP Solves the matrix completion problem using primal and dual
% SDP formulations with both IPNSDP and SeDuMi solvers.



clc; clear;

%% ==================== PRIMAL SDP FORMULATION ====================
% The primal problem minimizes the nuclear norm of a matrix X, which is
% a convex relaxation for the rank minimization problem.
% The problem is formulated as a Semi-Definite Program (SDP):
%
% min  0.5 * trace(Y)
% s.t. Y = [W, X; X', S] is Positive Semidefinite (PSD)
%      X(Omega) = M_known
%
% where Y is an (m+n)x(m+n) block matrix, and X is the m x n matrix to be completed.

%% ==================== DUAL SDP FORMULATION ====================
% The dual problem is a maximization problem over a dual variable matrix Z.
%
% max  sum_{i,j in Omega} (M_ij * Z_ij)
% s.t. [I, Z; Z', I] is Positive Semidefinite (PSD)
%      Z_ij = 0 for all (i,j) not in Omega
%
% where Z is the dual variable matrix of size m x n.

%% --------------------- Problem Data ---------------------
fprintf('Setting up matrix completion data...\n');
% Dimensions of the true low-rank matrix
m = 6;
n = 10;
% Generate a true low-rank matrix M_true
U = randn(m, 2);
V = randn(n, 2);
M_true = U * V';
% Create a mask of known entries (Omega)
p = 0.5;
Omega = find(rand(m, n) < p);
M_known = M_true(Omega);
% Create the incomplete matrix
M_incomplete = nan(m, n);
M_incomplete(Omega) = M_known;
fprintf('\nOriginal Matrix M_true:\n');
disp(M_true);
fprintf('\nIncomplete Matrix (NaN for missing):\n');
disp(M_incomplete);


%% ==================== PRIMAL SDP (IPNSDP) ====================
fprintf('\n------------------------------------------------\n');
fprintf('Solving Primal SDP with IPNSDP...\n');
fprintf('------------------------------------------------\n');
prob_p.name = 'mc_primal_sdp';
prob_p.problem_data = struct('m', m, 'n', n, 'M_known', M_known, 'Omega', Omega);
prob_p.nX = 1;
prob_p.dimX = [m + n];
prob_p.nu = 0;
prob_p.nPSDcon = 0;
prob_p.f_obj = @(x) 0.5 * trace(x.X{1});
% Use a local function to return the constraints in the required format {eq, ineq, psd}
prob_p.c1 = @(x) get_primal_constraints(x, prob_p);
prob_p.x0.X = {eye(m + n)};
prob_p.x0.u = [];
prob_p.solve = @() ipnsdp_solve(prob_p);
[xp_ipnsdp, ~] = prob_p.solve();
X_completed_ipnsdp = xp_ipnsdp.X{1}(1:m, m+1:end);
fprintf('\n--- IPNSDP Primal Results ---\n');
disp(X_completed_ipnsdp);
fprintf('Primal error (Frobenius norm): %.4f\n', norm(X_completed_ipnsdp - M_true, 'fro'));
fprintf('Rank: %d\n', rank(X_completed_ipnsdp, 1e-4));

%% ==================== DUAL SDP (IPNSDP) ====================
fprintf('\n------------------------------------------------\n');
fprintf('Solving Dual SDP with IPNSDP...\n');
fprintf('------------------------------------------------\n');
prob_d.name = 'mc_dual_sdp_modified';
prob_d.problem_data = struct('m', m, 'n', n, 'M_known', M_known, 'Omega', Omega);
prob_d.nX = 0;
prob_d.nu = m * n;
prob_d.dimX = [];
prob_d.nPSDcon = 1;
linear_idx_omega = Omega;
total_indices = 1:(m * n);
linear_idx_omega_complement = setdiff(total_indices, linear_idx_omega);
% Add linear_idx_omega_complement to problem data for local function access
prob_d.problem_data.linear_idx_omega_complement = linear_idx_omega_complement;
prob_d.f_obj = @(x) -sum(x.u(linear_idx_omega) .* M_known);
% Use a local function to return the constraints in the required format {eq, ineq, psd}
prob_d.c1 = @(x) get_dual_constraints(x, prob_d);
prob_d.x0.u = zeros(m * n, 1);
prob_d.solve = @() ipnsdp_solve(prob_d);
[xd_ipnsdp, ~] = prob_d.solve();
Z_ipnsdp = reshape(xd_ipnsdp.u, m, n);
fprintf('\n--- IPNSDP Dual Results ---\n');
disp(Z_ipnsdp);
fprintf('Optimal dual value: %.4f\n', -prob_d.f_obj(xd_ipnsdp));
fprintf('Check for sparsity: Max value outside Omega: %.4e\n', max(abs(Z_ipnsdp(linear_idx_omega_complement))));

%% ==================== PRIMAL SDP (SeDuMi) ====================
fprintf('\n\n------------------------------------------------\n');
fprintf('Solving Primal Problem with Sedumi...\n');
fprintf('------------------------------------------------\n');
fprintf('Generating Sedumi inputs...\n');
[At_primal, b_primal, c_primal, K_primal] = mc_primal_sedumi_full_vec(m, n, M_known, Omega);
fprintf('Calling Sedumi solver...\n');
[x_primal_sedumi, y_primal_sedumi, info_primal] = sedumi(At_primal, b_primal, c_primal, K_primal);
N_primal = m + n;
Y_primal_sedumi = reshape(x_primal_sedumi, N_primal, N_primal);
X_completed_sedumi = Y_primal_sedumi(1:m, m+1:end);
fprintf('\n--- Sedumi Primal Results ---\n');
fprintf('Completed Matrix X:\n');
disp(X_completed_sedumi);
primal_obj_sedumi = c_primal' * x_primal_sedumi;
fprintf('Optimal Primal Objective Value: %.4f\n', primal_obj_sedumi);
primal_error_fro = norm(X_completed_sedumi - M_true, 'fro');
fprintf('Primal completion error (Frobenius norm): %.4f\n', primal_error_fro);

%% ==================== DUAL SDP (SeDuMi) ====================
fprintf('\n------------------------------------------------\n');
fprintf('Solving Dual SDP with SeDuMi...\n');
fprintf('------------------------------------------------\n');
[At_d, b_d, c_d, K_d] = mc_dual_sedumi_full_vec(m, n, M_known, Omega);
[xd_sedumi, yd_sedumi, info_d] = sedumi(At_d, b_d, c_d, K_d);
Z_sedumi_vec = zeros(m * n, 1);
% The variable yd_sedumi only contains entries for Omega, so the
% constraint Z_ij=0 for (i,j) not in Omega is implicitly satisfied
% by the problem formulation.
Z_sedumi_vec(Omega) = yd_sedumi;
Z_sedumi = reshape(Z_sedumi_vec, m, n);
fprintf('\n--- Sedumi Dual Results ---\n');
disp(Z_sedumi);
dual_val = b_d' * yd_sedumi;
fprintf('Optimal dual value: %.4f\n', dual_val);
fprintf('Check for sparsity: Max value outside Omega: %.4e\n', max(abs(Z_sedumi(setdiff(1:(m*n), Omega)))));
end

%% ---------------- Local Functions ----------------
function [eq, ineq, psd] = get_primal_constraints(x, prob_p)
    % This function calculates the primal constraints for the IPNSDP solver.
    % It returns equality constraints, empty inequality constraints, and an empty PSD cell.
    eq = primal_eq_constraints(x, prob_p);
    ineq = [];
    psd = {}; % Empty cell array for PSD constraints
end

function [eq, ineq, psd] = get_dual_constraints(x, prob_d)
    % This function calculates the dual constraints for the IPNSDP solver.
    % It returns equality constraints, empty inequality constraints, and the PSD constraint in a cell.
    eq = dual_eq_constraints(x, prob_d);
    ineq = [];
    psd = {dual_psd_constraint_ipnsdp(x, prob_d)};
end

function c_eq = primal_eq_constraints(x, prob)
    % Enforces the primal equality constraints, Y(i, m+j) = M_known(k)
    % for all (i,j) in Omega.
    Y = x.X{1};
    m = prob.problem_data.m;
    n = prob.problem_data.n;
    Omega = prob.problem_data.Omega;
    M_known = prob.problem_data.M_known;
    [rows, cols] = ind2sub([m, n], Omega);
    num = length(Omega);
    c_eq = casadi.SX.zeros(num, 1);
    for k = 1:num
        c_eq(k) = Y(m + cols(k), rows(k)) - M_known(k);
    end
end

function Y = dual_psd_constraint_ipnsdp(x, prob)
    % Defines the Positive Semidefinite (PSD) cone constraint for the dual problem.
    % The block matrix [I, Z; Z', I] must be PSD.
    Z = reshape(x.u, prob.problem_data.m, prob.problem_data.n);
    Y = [casadi.SX.eye(prob.problem_data.m), Z; Z', casadi.SX.eye(prob.problem_data.n)];
end

function c_eq = dual_eq_constraints(x, prob)
    % Enforces the dual equality constraint that Z_ij = 0 for (i,j) not in Omega.
    linear_idx_omega_complement = prob.problem_data.linear_idx_omega_complement;
    c_eq = x.u(linear_idx_omega_complement);
end

function [At, b, c, K] = mc_dual_sedumi_full_vec(m, n, M_known, Omega)
    % MC_DUAL_SEDUMI_FULL_VEC generates SeDuMi inputs for the dual matrix completion problem.
    % The dual variable (y in SeDuMi) corresponds to the Z_ij values for (i,j) in Omega.
    % The constraint Z_ij=0 for (i,j) not in Omega is implicitly handled
    % by the size of the dual variable vector.
    num_vars = length(Omega);
    b = M_known(:);
    N_lmi = m + n;
    c = [eye(m), zeros(m, n); zeros(n, m), eye(n)];
    c = c(:);
    At = spalloc(num_vars, N_lmi^2, 2 * num_vars);
    [rows, cols] = ind2sub([m, n], Omega);
    for k = 1:num_vars
        i = rows(k); j = cols(k);
        idx_ur = (m + j - 1) * N_lmi + i;
        idx_ll = (i - 1) * N_lmi + (m + j);
        At(k, idx_ur) = 1;
        At(k, idx_ll) = 1;
    end
    K.s = N_lmi;
end

function [At, b, c, K] = mc_primal_sedumi_full_vec(m, n, M_known, Omega)
    % MC_PRIMAL_SEDUMI_FULL_VEC generates SeDuMi inputs for the primal matrix completion problem.
    % It formulates the problem as min c'x s.t. Ax=b, x in K, where x is the vectorized
    % block matrix Y and the constraints enforce the known entries.
    N = m + n;
    num_vars_primal = N * N;
    C_obj = 0.5 * eye(N);
    c = C_obj(:);
    b = M_known(:);
    num_cons = length(Omega);
    At = spalloc(num_vars_primal, num_cons, num_cons);
    [rows, cols] = ind2sub([m, n], Omega);
    for k = 1:num_cons
        i = rows(k);
        j = cols(k);
        primal_idx_ur = (m + j - 1) * N + i;
        At(primal_idx_ur, k) = 1;
    end
    K.s = N;
end