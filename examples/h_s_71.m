function h_s_71()
% H_S_71 solves the Hock and Shittkowski 71 without PSD
%   minimize    x1*x4*(x1 + x2 + x3) + x3
%   subject to
%       x1*x2*x3*x4 - 25 >= 0
%       x1^2 + x2^2 + x3^2 + x4^2 - 40 = 0
%       1 <= xi <= 5,   i = 1,...,4


clc; clear;

%% Problem Definition
prob.name = 'Hock_Shittkowski_71';
prob.problem_data = struct();

% Variables: no matrix X_k, all scalars in u
prob.nX = 0;
prob.dimX = [];
prob.nu = 4;   % u = [x1; x2; x3; x4]
prob.lbX = [];
prob.lbu = ones(4,1);
prob.ubu = 5*ones(4,1);

% nPSD constraints
prob.nPSDcon = 0;       

% Objective
prob.f_obj = @(x) objective(x, prob);

% Nonlinear and PSD constraints
prob.c1 = @(x) nonlinear_constraint(x, prob);

% Initial guess
prob.x0.X = {};
prob.x0.u = 1*[2; 2; 2; 2]; % Feasible start

%% Solve
prob.solve = @() ipnsdp_solve(prob);
[x_sol, info] = prob.solve();

disp('Solution:');
disp('x ='), disp(x_sol.u);
disp('Objective value:'), disp(info.obj_value);

end

%% Objective Function
function f_val = objective(x, prob)
u = x.u; 
f_val = u(1)*u(4)*(u(1) + u(2) + u(3)) + u(3);
end

%% Nonlinear & PSD Constraints
function [eq, ineq, psd] = nonlinear_constraint(x, prob)
u = x.u; % [x1; x2; x3; x4]

% eq = []; % Equalities

eq = u(1)^2 + u(2)^2 + u(3)^2 + u(4)^2 - 40;

% Inequalities (<= 0)
ineq = - [
    u(1)*u(2)*u(3)*u(4) - 25;
];

% No PSD matrix
psd = {};

end
