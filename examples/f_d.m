function f_d()
% F_D solves a financial discrimination problem using a Logit model
% with a positive semidefinite quadratic discriminant function.
% The problem is formulated for the IPNSDP solver.
% The maximization objective is converted to a minimization problem by negation.
%
% Problem:
%   maximize_{a, b, Q}  sum_{i=1 to M} (y_i * z(x_i) - log(1 + e^{z(x_i)}))
%   subject to: Q >= 0
%   where z(x) = a + bᵗx + 0.5 * x^TQx
%
% The solver assumes the matrix variable (X) is in the PSD cone.

clc; clear;

q_values = [6, 8, 10, 12 ];
save_folder = 'sol';

% Create folder if it doesn't exist
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for i = 1:length(q_values)
    %% Problem Data

    q = q_values(i);

    %% --------------------- Problem Data ---------------------
    M = 12000;       % Number of companies (data points)
    % Define parameter ranges
    mu_min = 10;
    mu_max = 1000;
    sigma_min = 1;
    sigma_max = 200;

    % Create random mu and sigma
    mu = mu_min + (mu_max - mu_min) .* rand(1, q);
    sigma = sigma_min + (sigma_max - sigma_min) .* rand(1, q);

    % Generate raw normal data
    x_raw = randn(M, q);

    % Scale and shift data
    x_data = x_raw .* repmat(sigma, M, 1) + repmat(mu, M, 1);

    % Ensure non-negative values
    x_data = max(0, x_data)/(max(max(x_data)));
    % x_data = randn(M, q); 
    
    % Generate random binary labels (0: non-failure, 1: failure)
    y_data = randi([0, 1], M, 1);
    
    %% Problem Setup 
    % Problem structure
    prob.name = 'financial_discrimination';
    prob.problem_data = struct('M', M, 'q', q, 'x_data', x_data, 'y_data', y_data);
    
    % Variables
    prob.nX = 1;               % One matrix variable: Q
    prob.dimX = [q];           % Dimension of matrix Q
    prob.nu = q + 1;           % Vector variables: b (q) and a (1)

    % nPSD constraints
    prob.nPSDcon = 0;          %  PSD constraint for Q is handled by the X block

    % Objective function handle
    prob.obj = @(x) get_objective_value(x, prob);
    
    % Constraint function handle
    prob.nlcon = @(x) get_constraints_financial(x, prob);
    
    % Initial guesses
    prob.x0.X = {1* eye(q)};           % Initial Q
    prob.x0.u = .01*ones(q + 1, 1);     % Initial [a; b]

    % Solve the problem using the IPNSDP solver
    [x, info] = ipnsdp_solve(prob);
    
    % Extract solution
    Q_sol = x.X{1};        % Optimal Q
    a_sol = x.u(1);        % Optimal a
    b_sol = x.u(2:end);    % Optimal b
    
    % Example: Evaluate discriminant on a new point x_new
    % x_new = randn(q, 1);
    % z_new = a_sol + b_sol' * x_new + 0.5 * x_new' * Q_sol * x_new;
    % p = 1 / (1 + exp(-z_new));  % Probability of failure

    %% Save result

    filename = sprintf('%s/f_d_q_%d.mat', save_folder, q);
    save(filename, "prob", "x", "info");

end
end

function obj_val = get_objective_value(x_solver_vars, prob)
    % Computes the negative Logit objective value for minimization.
    y_data = prob.problem_data.y_data;
    z_x = get_z_x(x_solver_vars, prob);
    obj_val = -sum(y_data .* z_x - log(1 + exp(z_x)));
end

function z_x = get_z_x(x, prob)
    % Computes z(x_i) = a + bᵗx_i + 0.5 * x_iᵗQx_i for all data points.
    q = prob.problem_data.q;
    x_data = prob.problem_data.x_data;
    
    Q = x.X{1};
    a = x.u(1);
    b = x.u(2:end);

    z_x = casadi.SX.zeros(size(x_data, 1), 1);
    for i = 1:size(x_data, 1)
        xi = x_data(i, :)';
        z_x(i) = a + b' * xi + 0.5 * xi' * Q * xi;
    end
end

function [eq, ineq, psd] = get_constraints_financial(~, ~)
    % No additional equality or inequality constraints in this model.
    eq = [];
    ineq = [];
    psd = {}; % PSD handled directly by Q \in S^+ in the X block
end
