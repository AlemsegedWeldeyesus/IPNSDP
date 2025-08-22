function bmi_1()
% BMI_1 solves a Bilinear Matrix Inequality (BMI) problem using IPNSDP.
%
% Problem formulation:
%   minimize_{x, lambda}  lambda
%   subject to:
%       -c_l <= x_l <= c_l,  for l = 1,...,n
%       B(x) = A0 + sum_{k=1}^n x_k*A_k + sum_{i=1}^d sum_{j=1}^d x_i*x_j*K_{i,j} >= lambda * I_m
%
% A_k and K_{i,j} are sparse symmetric matrices generated using sprandn.

clc; clear;

ndm_values = [50 25 50; 50, 50, 50; 100, 50, 100; 100, 100, 100];
save_folder = 'sol';

% Create folder if it doesn't exist
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for i =1:size(ndm_values,1)
    n = ndm_values(i,1); % Dimension of x
    d = ndm_values(i,2); % Dimension of quadratic terms (d > n^2)
    m = ndm_values(i,3); % Dimension of symmetric matrices A_k and K_{i,j}
    %% --------------------- Problem Data ---------------------
    r = 0.2;  % Sparsity factor for random matrices
    
    % Generate box constraint bounds
    c = abs(randn(n, 1));
    
    % Generate sparse symmetric matrix A0
    A0 = symmetrize(sprandn(m, m, r));
    
    % Generate A_k matrices for linear terms
    A_k = cell(n, 1);
    for k = 1:n
        A_k{k} = symmetrize(sprandn(m, m, r));
    end
    
    % Generate K_{i,j} matrices for bilinear terms
    K_ij = cell(d, d);
    for i = 1:d
        for j = 1:d
            K_ij{i,j} = symmetrize(sprandn(m, m, r));
        end
    end
    
    %% Problem Setup 
    prob.name = 'bmi_primal_sdp';
    
    % Store problem constants
    prob.problem_data = struct( ...
        'n', n, ...
        'd', d, ...
        'm', m, ...
        'c_abs', c, ...
        'A0', A0, ...
        'A_k', {A_k}, ...
        'K_ij', {K_ij} ...
    );
    
    % Variables
    prob.nX = 0;             % No matrix variables
    prob.dimX = [];          
    prob.nu = n + 1;         % Variables: x (n) and lambda (1)
    prob.lbu= [-c;-inf];
    prob.ubu= [c;inf];

    % nPSD constraints
    prob.nPSDcon = 1;        % One PSD constraint (the BMI)
    
    % Objective: minimize lambda (last variable)
    prob.f_obj = @(x) x.u(end);

    % Initial guess
    prob.x0.u = zeros(n + 1, 1);
    
    % Nonlinear constraints
    prob.c1 = @(x) get_constraints_bmi(x, prob);
    
    % Solve 
    prob.solve = @() ipnsdp_solve(prob);
    [x, info] = prob.solve();

    %% Save result
    filename = sprintf('%s/bmi_n%d_d%d_m%d.mat', save_folder, n,d,m);
    save(filename, "prob", "x", "info");

end
end

function [eq, ineq, psd] = get_constraints_bmi(x, prob)
% Returns constraint evaluations: equality, inequality, and PSD constraints

    u = x.u(1:end-1);       % Decision variables
    lambda = x.u(end);      % Lambda variable

    n = prob.problem_data.n;
    d = prob.problem_data.d;
    m = prob.problem_data.m;
    c = prob.problem_data.c_abs;
    A0 = prob.problem_data.A0;
    A_k = prob.problem_data.A_k;
    K_ij = prob.problem_data.K_ij;

    % No equality constraints
    eq = [];

    % No inequality constraints
    ineq =[];

    % PSD constraint: B(x) - lambda*I_m >= 0
    Bx = A0;
    for k = 1:n
        Bx = Bx + u(k) * A_k{k};
    end
    for i = 1:d
        for j = 1:d
            Bx = Bx + u(i) * u(j) * K_ij{i,j};
        end
    end

    psd = { -Bx + lambda * eye(m) };
end

function S = symmetrize(M)
% Ensures matrix M is symmetric
    S = (M + M') / 2;
end
