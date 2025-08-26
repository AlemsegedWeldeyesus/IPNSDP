function m_m_ev()
% M_M_EV solves the minimization of the minimum eigenvalue problem
% using the IPNSDP solver.
%
% Problem:
%   minimize     trace(Pi * M(q))
%   subject to   trace(Pi) = 1
%                Pi >= 0
%                -1 <= q_i <= 1  for i = 1, 2, 3
%
% where M(q) = q1*M1 + q2*M2 + q3*M3.
% M1, M2, M3 are symmetric matrices with random entries in [-1, 1].

clc; clear;

%% --------------------- Problem Data ---------------------
n = 110; % Dimension of the symmetric matrices Pi and M(q)
p = 3; % Number of components for the vector q
% Generate random symmetric matrices M1, M2, M3
M1_rand = rand(n) * 2 - 1;
M1 = (M1_rand + M1_rand') / 2;
M2_rand = rand(n) * 2 - 1;
M2 = (M2_rand + M2_rand') / 2;
M3_rand = rand(n) * 2 - 1;
M3 = (M3_rand + M3_rand') / 2;

%% ==================== IPNSDP FORMULATION ====================
prob_p.name = 'min_eigenvaluenlcon_primal_sdp';
% Store problem constants in the problem_data struct
prob_p.problem_data = struct('n', n, 'p', p, 'M1', M1, 'M2', M2, 'M3', M3);
prob_p.nX = 1;      % One matrix variable: Pi
prob_p.dimX = [n];  % The dimension of the matrix variable
prob_p.nu = p;      % Three vector variables: q1, q2, q3
prob_p.lbu = -ones(p,1);
prob_p.ubu = ones(p,1);

% nPSD constraints
prob_p.nPSDcon = 0; % PSD constraint for pi is handled by the X block

% Objective function
prob_p.obj = @(x) trace(x.X{1} * get_M_q(x, prob_p));

% Solve the problem
prob_p.nlcon = @(x) get_constraints_min_eigenvalue(x, prob_p);

% Initial guess for the variables
prob_p.x0.X = {eye(n)}; % Initial guess for Pi
prob_p.x0.u = zeros(p, 1); % Initial guess for q

% Solve
[x_sol, info] = ipnsdp_solve(prob_p);

end

function [eq, ineq, psd] = get_constraints_min_eigenvalue(x_solver_vars, prob_p)
    % Extract variables from the solver's single vector and matrix variable
    n = prob_p.problem_data.n;
    q_vars = x_solver_vars.u;
    Pi_var = x_solver_vars.X{1};

    % Equality constraint: trace(Pi) = 1
    eq = trace(Pi_var) - 1;

    % No inequality 
    ineq = [];

    % No extra PSD constraints (Pi >= 0 is handled in the X block)
    psd = {};
end

function M_q = get_M_q(x, prob_p)
    % This function calculates the matrix M(q)
    M1 = prob_p.problem_data.M1;
    M2 = prob_p.problem_data.M2;
    M3 = prob_p.problem_data.M3;
    q = x.u;
    q1 = q(1);
    q2 = q(2);
    
    M_q = q1 * q2* M1 + q1 * M2 + q2 * M3;
end
