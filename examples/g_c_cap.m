function g_c_cap()
% G_C_CAP solves the Gaussian channel capacity problem using IPNSDP.
% The objective is:
%
%   maximize_{x, t} (1/2) * sum_{j=1 to N} log(1 + t_j)
%
% subject to:
%   (1/N) * sum_{j=1 to N} x_j <= 1
%   x_j >= 0, t_j >= 0
%   [ 1 - a_j * t_j, sqrt(r_j);
%     sqrt(r_j),     a_j * x_j + r_j ] >= 0, for j = 1,...,N
%
% Note:
% - Reformulated as a minimization problem for IPNSDP.
% - The variables x and t are constrained to be non-negative, so |.| is removed.

clc; clear;
clc; clear;

N_values = [100, 200, 400, 800, 1600, 3200 ];
save_folder = 'sol';

% Create folder if it doesn't exist
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for i = 1:length(N_values)
    % Problem size
    N = N_values(i);

    % Random channel constants a_j and noise powers r_j in [0, 1]
    a_j = rand(N, 1);
    r_j = rand(N, 1);
    
    % Define problem structure
    prob.name = 'gaussian_channel_capacity_primal_sdp';
    prob.problem_data = struct('N', N, 'a_j', a_j, 'r_j', r_j);
    
    % Variables
    prob.nX = 0;            % No matrix variables
    prob.nu = 2 * N;        % Variables: x_1,...,x_N and t_1,...,t_N
    prob.dimX = [];
    prob.lbu=zeros(2*N,1);
    
    % nPSD constraints
    prob.nPSDcon = N;       % N: One 2x2 PSD constraint per channel
    
    % Objective function (negated for minimization)
    % x.u = [x_1,...,x_N, t_1,...,t_N]
    prob.f_obj = @(x) -0.5 * sum(log(1 + x.u(N+1:end)));
    
    % Constraints: equality, inequality, and PSD
    prob.c1 = @(x) get_constraints_channel_capacity(x, prob);
    
    % Initial guess: concatenate x and t as a 2N-dimensional vector
    prob.x0.u = ones(2 * N, 1);

    % Call the IPNSDP solver
    prob.solve = @() ipnsdp_solve(prob);
    [x, info] = prob.solve();
    
     %% Save result

    filename = sprintf('%s/g_c_cap_N%d.mat', save_folder, N);
    save(filename, "prob", "x", "info");
end




end

function [eq, ineq, psd] = get_constraints_channel_capacity(x_solver_vars, prob)
% GET_CONSTRAINTS_CHANNEL_CAPACITY generates constraints for the Gaussian channel capacity problem.
%
% Inputs:
%   - x_solver_vars: structure with field u (concatenated x and t variables)
%   - prob: problem structure containing a_j, r_j, N
%
% Outputs:
%   - eq: equality constraints (empty)
%   - ineq: inequality constraints (vector)
%   - psd: cell array of PSD constraints (length N, each 2x2)

    % Extract problem data and variables
    N = prob.problem_data.N;
    a_j = prob.problem_data.a_j;
    r_j = prob.problem_data.r_j;

    u = x_solver_vars.u;
    x_vars = u(1:N);
    t_vars = u(N+1:end);

    % No equality constraints
    eq = [];

    % Inequality constraints (IPNSDP expects form <= 0):
    ineq = [(1/N) * sum(x_vars) - 1];

    % PSD constraints:
    % For each j, enforce:
    % [ 1 - a_j * t_j, sqrt(r_j);
    %   sqrt(r_j),     a_j * x_j + r_j ] >= 0
    psd = cell(N, 1);
    for j = 1:N
        psd{j} = [1 - a_j(j) * t_vars(j), sqrt(r_j(j));
                  sqrt(r_j(j)),           a_j(j) * x_vars(j) + r_j(j)];
    end

end
