function h_s_43()
% H_S_43 solves the Hock and Shittkowski 43 (Rosen-Suzuki Problem) without PSD
%   minimize    x1^2 + x2^2 + 2*x3^2 + x4^2 - 5*x1 - 5*x2 - 21*x3 + 7*x4
%   subject to
%       8 - x1^2 - x2^2 - x3^2 - x4^2 - x1 + x2 - x3 + x4 >= 0
%       10 - x1^2 - 2*x2^2 - x3^2 - 2*x4^2 + x1 + x4 >= 0
%       5 - 2*x1^2 - x2^2 - x3^2 - 2*x1 + x2 + x4 >= 0


%% Problem Definition
prob.name = 'Hock_Shittkowski_43';
prob.problem_data = struct();

% Variables: no matrix X_k, all scalars in u
prob.nX = 0;
prob.dimX = [];
prob.nu = 4;   % u = [x1; x2; x3; x4]
prob.lbX = [];

% nPSD constraints
prob.nPSDcon = 0;        

% Objective
prob.obj = @(x) objective(x, prob);

% Nonlinear and PSD constraints
prob.nlcon = @(x) nonlinear_constraint(x, prob);

% Initial guess
prob.x0.X = {};
prob.x0.u = 10*[1; 1; 1; 1];

%% Solve
[x_sol, info] = ipnsdp_solve(prob);

disp('Solution:');
disp('x ='), disp(x_sol.u);
disp('Objective value:'), disp(info.obj_value);

end

%% Objective Function
function f_val = objective(x, prob)
u = x.u; % [x1; x2; x3; x4]
f_val = u(1)^2 + u(2)^2 + 2*u(3)^2 + u(4)^2 ...
        - 5*u(1) - 5*u(2) - 21*u(3) + 7*u(4);
end

%% Nonlinear & PSD Constraints
function [eq, ineq, psd] = nonlinear_constraint(x, prob)
u = x.u; % [x1; x2; x3; x4]

% Inequalities (<= 0)
ineq = - [
    8  - u(1)^2 - u(2)^2 - u(3)^2 - u(4)^2 - u(1) + u(2) - u(3) + u(4);
    10 - u(1)^2 - 2*u(2)^2 - u(3)^2 - 2*u(4)^2 + u(1) + u(4);
    5  - 2*u(1)^2 - u(2)^2 - u(3)^2 - 2*u(1) + u(2) + u(4)
];

eq = []; % No equalities

% No PSD matrix
psd = {}; 

end
