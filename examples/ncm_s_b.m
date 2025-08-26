function ncm_s_b()
%NCM_S_B Solves the scaled Nearest Correlation Matrix (NCM) problem:
%
%   minimize    || z * X - H ||_F^2
%   subject to  diag(z * X) = 1
%               I <= X <= κI
%               X symmetric
%
% This problem is formulated and solved using the ipnsdp solver.

clc; clear;

Hsize = [10, 20, 40, 80];
save_folder = 'sol';

% Create folder if it doesn't exist
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for i = 1:length(Hsize)
    %% Problem Data

    n = Hsize(i);

    % Generate a random correlation matrix X_true
    A = randn(n);
    X_true = A * A' + 0.2 * eye(n);

    % Add symmetric noise
    E = randn(n);
    E = (E + E') / 2;
    H = X_true + 0.2 * E;

    % Set diagonal to 1
    H(1:n+1:end) = 1;
    
    kappa = 10;

    %% Problem Definition

    prob.name = 'ncm_s';
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
    prob.obj = @(x) norm(x.u * x.X{1} - H, 'fro')^2;
  
    % Constraints: [eq, ineq, psd]
    prob.nlcon = @(x) deal( ...
        diag(x.u * x.X{1}) - 1, ...       % diag(z*X) = 1
        [], ...                           % no inequalities
        {});                              % no PSD constraints
    
    % Initial point
    % prob.x0.X = { 0.5 * kappa * eye(n) }; % Feasible start
    % prob.x0.u = 1;

    %% Solve
    [x, info] = ipnsdp_solve(prob);

    %% Save result
    filename = sprintf('%s/ncm_s_b_%d.mat', save_folder, n);
    save(filename, "prob", "x", "info");
end

end
