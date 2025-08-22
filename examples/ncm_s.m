function ncm_s()
%NCM_S Solves the scaled Nearest Correlation Matrix (NCM) problem:
%
%   minimize    || z * X - H ||_F^2
%   subject to  diag(z * X) = 1
%               I <= X <= κI
%               X symmetric
%
% This problem is formulated and solved using the ipnsdp solver.

clc; clear;

%% Problem Data
H = [  1.00, -0.44, -0.20,  0.81, -0.46, -0.05;
      -0.44,  1.00,  0.87, -0.38,  0.81, -0.58;
      -0.20,  0.87,  1.00, -0.17,  0.65, -0.56;
       0.81, -0.38, -0.17,  1.00, -0.37, -0.15;
      -0.46,  0.81,  0.65, -0.37,  1.00,  0.08;
      -0.05, -0.58, -0.56, -0.15,  0.08,  1.00 ];

kappa = 10;
n = size(H, 1);

%% Problem Definition
prob.name         = 'ncm_s';
prob.problem_data = struct('H', H, 'kappa', kappa);

% Problem dimensions
prob.nX      = 1;                     % One matrix variable X
prob.dimX    = [n];                   % X in S^n
prob.nu      = 1;                     % One scalar variable z
prob.lbX     = [1];                   % Enforce X >= I
prob.ubX     = [kappa];               % Enforce X <= kappa*I

% nPSD constraints
prob.nPSDcon = 0;        

% Objective function
prob.f_obj = @(x) norm(x.u * x.X{1} - H, 'fro')^2;

% Constraints: [eq, ineq, psd]
prob.c1 = @(x) deal( ...
    diag(x.u * x.X{1}) - 1, ...       % diag(z*X) = 1
    [], ...                           % no inequalities
    {});                              % no PSD constraints

% Initial point
% prob.x0.X = { 0.5 * kappa * eye(n) }; % Feasible start
% prob.x0.u = 1;

%% Solve
prob.solve = @() ipnsdp_solve(prob);
[x, info] = prob.solve();

%% Results
fprintf('Optimal z*X:\n');
disp(x.u * x.X{1});

fprintf('Optimal z^{-1} = %.6f\n', 1 / x.u);
fprintf('Final objective value: %.6f\n', info.obj_value);

end
