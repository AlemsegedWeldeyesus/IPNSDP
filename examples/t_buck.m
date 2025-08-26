function t_buck()
% This is a minimization problem from the paper "Truss topology optimization with global stability
% constraints" by Weldeyesus et al. The objective is to minimize l' * a
% subject to the following constraints:
% 1. f' * u <= zeta
% 2. K(a) * u = f
% 3. K(a) + tau * G(a, u) >= 0
% 4. a >= 0
% 5. u is a free variable
% Note that v is fixed, therefore no geometry optimization.
clc
clear

% Import problem data
problem_data = getProblemData();

% Define the IPNSDP problem structure
prob = struct();
prob.name = 'GeoOpt_Buckling';
prob.nX = 0;
prob.dimX = [];
prob.nu = problem_data.n;
prob.nPSDcon = problem_data.nL;
prob.lbu = problem_data.lbu;
prob.A21 = problem_data.A21;
prob.A22 = problem_data.A22;
prob.b2  = problem_data.b2;
prob.x0.X = cell(prob.nX, 1);
prob.x0.u = problem_data.u0;
prob.problem_data = problem_data;
prob.obj = @(x) geopt_objective(x, prob);
prob.nlcon = @(x) geopt_nonlinear_con(x, prob);

% Solve
[x_sol, info] = ipnsdp_solve(prob);

% Visualize
ipmsgeopt_showTruss(x_sol, prob);

end


function problem_data = getProblemData()

problem_data.dim = 3;
problem_data.box = [8, 1, 2];
problem_data.divx = 16/2*1;
problem_data.divy = 2/2*1;
problem_data.divz = 2*1;

problem_data.force = [
    1,0,0,0,0,-.001, 1;
    2,0,0,0,0,-.001, 1;
    3,0,0,0,0,-.001, 1;
    4,0,0,0,0,-.001, 1;
    5,0,0,0,0,-.001, 1;
    6,0,0,0,0,-.001, 1;
    7,0,0,0,0,-.001, 1;
    1,1,0,0,0,-.001, 1;
    2,1,0,0,0,-.001, 1;
    3,1,0,0,0,-.001, 1;
    4,1,0,0,0,-.001, 1;
    5,1,0,0,0,-.001, 1;
    6,1,0,0,0,-.001, 1;
    7,1,0,0,0,-.001, 1
];

problem_data.support = [
    0,0,0, 1,1,1;
    0,1,0, 1,1,1;
    8,0,0, 1,1,1;
    8,1,0, 1,1,1
];

problem_data.c = .003;
problem_data.r = 0.2;
problem_data.LoadingFcator = ones(max(problem_data.force(:,end)), 1);

[problem_data.nL, problem_data.Nd, problem_data.NdFixed] = getNd(problem_data);
[problem_data.dof, problem_data.dofgeo, problem_data.nodegeo] = getDof(problem_data);
[problem_data.Cn, problem_data.Cnlst] = getCnCnlst(problem_data);
problem_data.f = getForce(problem_data);

nDispDofs = nnz(problem_data.dof) / problem_data.nL;
nTotalDispDofs = nnz(problem_data.dof);
nGeoDofs = nnz(problem_data.dofgeo);
nMembers = length(problem_data.Cnlst);
nL = problem_data.nL;

problem_data.A21 = zeros(nL + nMembers, 0);

f = zeros(nL, size(problem_data.f,1));
nf = size(problem_data.f,1) / nL;
for i = 1:nL
    f(i, (i-1)*nf+1 : i*nf) = problem_data.f((i-1)*nf+1 : i*nf);
end

problem_data.A22 = [
    sparse(nL, nGeoDofs), zeros(nL, nMembers), f;
];

xyz0 = [problem_data.Nd(:,1:problem_data.dim)]';
xyz0 = xyz0(:);
xyz0(problem_data.dofgeo == 0) = [];

problem_data.b2 = [problem_data.c];

problem_data.X0 = cell(problem_data.nL,1);
for j = 1:problem_data.nL
    problem_data.X0{j} = 1e0 * eye(nTotalDispDofs);
end

problem_data.u0 = [xyz0; 10*ones(nMembers,1); zeros(problem_data.nL * nDispDofs,1)];
problem_data.m = problem_data.nL;
problem_data.n = nGeoDofs + nMembers + problem_data.nL * nDispDofs;
problem_data.lbu=[-inf*ones(nGeoDofs,1);zeros(nMembers,1); -inf*ones(problem_data.nL * nDispDofs,1)];
end


function [nL,Nd,NdFixed] = getNd(p)

dim = p.dim;
divx = p.divx;
divy = p.divy;
if dim == 3
    divz = p.divz;
end

force = p.force;
box = p.box;
support = p.support;
nL = max(force(:,end));

if dim == 3
    [x,y,z] = meshgrid(linspace(0,box(1),divx+1), linspace(0,box(2),divy+1), linspace(0,box(3),divz+1));
    Nd(:,1:3) = [x(:),y(:),z(:)];
else
    [x,y] = meshgrid(linspace(0,box(1),divx+1), linspace(0,box(2),divy+1));
    Nd(:,1:2) = [x(:),y(:)];
end

NdFixed = 0 * Nd(:,1:dim);

s_idx = findNode(Nd(:,1:dim), support(:,1:dim), 1e-6);
Nd(s_idx,dim+1:2*dim) = support(:,dim+1:2*dim);
NdFixed(s_idx,1:dim) = Nd(s_idx,1:dim);

for i = 1:nL
    f_idx_all = find(force(:,end) == i);
    f_idx = findNode(Nd(:,1:dim), force(f_idx_all,1:dim), 1e-6);
    Nd(f_idx,2*dim+(i-1)*dim+1 : 2*dim+i*dim) = force(f_idx_all,dim+1:2*dim);
    Nd(f_idx,2*dim+nL*dim+i) = force(f_idx_all,2*dim+1);
    NdFixed(f_idx,1:dim) = Nd(f_idx,1:dim);
end

p.MoreFixedNodes = Nd(:,1:dim);
if (isfield(p,'MoreFixedNodes'))
    mf_idx = findNode(Nd(:,1:dim), p.MoreFixedNodes,1e-6);
    NdFixed(mf_idx,1:dim)= p.MoreFixedNodes;
    Nd(mf_idx,2*dim+nL*dim+nL+1) = 1;
end

end


function [dof, dofgeo, nodegeo] = getDof(p)

Nd = p.Nd;
nL = p.nL;
dim = p.dim;
m = size(Nd, 1);

dof = ones(m*dim, 1);

if dim == 2
    dof(find(Nd(:,3) == 1)*2-1) = 0;
    dof(find(Nd(:,4) == 1)*2) = 0;
elseif dim == 3
    dof(find(Nd(:,4) == 1)*3-2) = 0;
    dof(find(Nd(:,5) == 1)*3-1) = 0;
    dof(find(Nd(:,6) == 1)*3) = 0;
end

dofgeo = dof;
dof = repmat(dof, nL, 1);

m2 = size(Nd, 2);
for j = 1:nL
    if dim == 2
        for i = 1:m
            if Nd(i, (j-1)*dim+5) ~= 0 || Nd(i, (j-1)*dim+6) ~= 0 || Nd(i, m2) ~= 0
                dofgeo(i*2-1) = 0;
                dofgeo(i*2) = 0;
            end
        end
    elseif dim == 3
        for i = 1:m
            if Nd(i, (j-1)*dim+7) ~= 0 || Nd(i, (j-1)*dim+8) ~= 0 || Nd(i, (j-1)*dim+9) ~= 0 || Nd(i, m2) ~= 0
                dofgeo(i*3-2) = 0;
                dofgeo(i*3-1) = 0;
                dofgeo(i*3) = 0;
            end
        end
    end
end

nodegeo = reshape(dofgeo, dim, [])';
nodegeo = nodegeo(:,1);

end


function [cn, cnlst] = getCnCnlst(p)

n = p.Nd;
dim = p.dim;

cn = zeros(size(n,1)*(size(n,1)-1)/2, 3);
count = 0;

for i = 1:size(n,1)-1
    for j = i+1:size(n,1)
        count = count + 1;
        l = norm(n(i,1:dim) - n(j,1:dim));
        cn(count,:) = [i, j, l];
    end
end

minstep = (max(n(:,1)) - min(n(:,1))) / p.divx;
cnlst = find(cn(:,3) < 1*sqrt(dim + 1e-3)*minstep);

end


function F = getForce(p)

Nd = p.Nd;
nL = p.nL;
dim = p.dim;
dof = p.dof;
m = size(Nd, 1);

F = zeros(m*dim*nL, 1);

for i = 1:m
    for j = 1:nL
        if (Nd(i, 2*dim + dim*nL + j) == j)
            F((i-1)*dim + (j-1)*m*dim + 1 : i*dim + (j-1)*m*dim) = ...
                Nd(i, 2*dim + (j-1)*dim + 1 : 2*dim + j*dim)';
        end
    end
end

F(dof == 0, :) = [];

% Return as column vector
F = F;

end


function obj = geopt_objective(x, prob)

p = prob.problem_data;
import casadi.*

mc = nnz(p.dofgeo);
n = length(p.Cnlst);
a = x.u(mc+1:mc+n);
Nd = getNewNd(p, x);
cn = p.Cn;
cnlst = p.Cnlst;

L = SX.zeros(size(a));

for i = 1:length(cnlst)
    cid = cnlst(i);
    n1 = cn(cid,1);
    n2 = cn(cid,2);
    l = norm(Nd(n2, 1:p.dim) - Nd(n1, 1:p.dim));
    L(i) = l;
end

obj = L' * a;

end


function [Ceq, Cinq, c_mat] = geopt_nonlinear_con(x, prob)

problem_data = prob.problem_data;
import casadi.*

Cinq = [];
m1 = nnz(problem_data.dof);
mc = nnz(problem_data.dofgeo);
n = length(problem_data.Cnlst);

a = x.u(mc+1:mc+n);
d1 = x.u(mc+n+1:mc+n+m1);
nL = problem_data.nL;
dof = problem_data.dof(1:length(problem_data.dof)/nL);

Nd = getNewNd(problem_data,x);
dim = problem_data.dim;
Cn = problem_data.Cn;
Cnlst = problem_data.Cnlst;

L = SX.zeros(size(a));

rows_B_full = size(Nd, 1) * dim;
cols_B = length(Cnlst);

B_full = SX.zeros(rows_B_full, cols_B);
Bd_full = SX.zeros(rows_B_full, cols_B);
Be_full = SX.zeros(rows_B_full, cols_B);

if dim == 2
    for i = 1:length(Cnlst)
        cid = Cnlst(i);
        n1 = Cn(cid, 1);
        n2 = Cn(cid, 2);
        xc = Nd(n2, 1) - Nd(n1, 1);
        y = Nd(n2, 2) - Nd(n1, 2);
        l = sqrt(xc^2 + y^2);
        c = xc / l;
        s = y / l;
        L(i) = l;
        
        B_full((n1-1)*2+1, i) = -c;
        B_full((n1-1)*2+2, i) = -s;
        B_full((n2-1)*2+1, i) = c;
        B_full((n2-1)*2+2, i) = s;
        
        Bd_full((n1-1)*2+1, i) = s;
        Bd_full((n1-1)*2+2, i) = -c;
        Bd_full((n2-1)*2+1, i) = -s;
        Bd_full((n2-1)*2+2, i) = c;
    end
elseif dim == 3
    for i = 1:length(Cnlst)
        cid = Cnlst(i);
        n1 = Cn(cid, 1);
        n2 = Cn(cid, 2);
        xc = Nd(n2, 1) - Nd(n1, 1);
        y = Nd(n2, 2) - Nd(n1, 2);
        z = Nd(n2, 3) - Nd(n1, 3);
        l = sqrt(xc^2 + y^2 + z^2);
        c = xc / l;
        s = y / l;
        cs = z / l;
        L(i) = l;
        
        B_full((n1-1)*3+1, i) = -c;
        B_full((n1-1)*3+2, i) = -s;
        B_full((n1-1)*3+3, i) = -cs;
        B_full((n2-1)*3+1, i) = c;
        B_full((n2-1)*3+2, i) = s;
        B_full((n2-1)*3+3, i) = cs;
        
        w = cs;
        [dx, dz] = my_direction_function(c, s, w);
        d11 = dx(1); d2 = dx(2); d3 = dx(3);
        e1 = dz(1); e2 = dz(2); e3 = dz(3);
        
        Bd_full((n1-1)*3+1, i) = -d11;
        Bd_full((n1-1)*3+2, i) = -d2;
        Bd_full((n1-1)*3+3, i) = -d3;
        Bd_full((n2-1)*3+1, i) = d11;
        Bd_full((n2-1)*3+2, i) = d2;
        Bd_full((n2-1)*3+3, i) = d3;
        
        Be_full((n1-1)*3+1, i) = -e1;
        Be_full((n1-1)*3+2, i) = -e2;
        Be_full((n1-1)*3+3, i) = -e3;
        Be_full((n2-1)*3+1, i) = e1;
        Be_full((n2-1)*3+2, i) = e2;
        Be_full((n2-1)*3+3, i) = e3;
    end
end

rows_to_keep = find(dof ~= 0);
B = B_full(rows_to_keep, :);
Bd = Bd_full(rows_to_keep, :);
Be = Be_full(rows_to_keep, :);

Kau = cell(nL, 1);
ak_tqG = cell(nL, 1);

for j = 1:nL
    dj = d1((j-1)*m1/nL + 1 : j*m1/nL);
    Kau{j} = B * diag(a ./ L) * B' * dj;
    if dim == 2
        ak_tqG{j} = B * diag(a ./ L) * B' + problem_data.LoadingFcator(j) * Bd * diag((a ./ (L.^2)) .* (B' * dj)) * Bd';
    elseif dim == 3
        ak_tqG{j} = B * diag(a ./ L) * B' + problem_data.LoadingFcator(j) * (Bd * diag((a ./ (L.^2)) .* (B' * dj)) * Bd' + Be * diag((a ./ (L.^2)) .* (B' * dj)) * Be');
    end
end

Kau2 = [];
for j = 1:nL
    Kj = Kau{j};
    Kau2 = [Kau2; Kj(:)];
end

Ceq = Kau2 - problem_data.f;

C2 = sum((Nd(:,1:dim) - problem_data.Nd(:,1:dim)).^2, 2);
nodegeo_indices_to_keep = find(problem_data.nodegeo == 1);
C2 = C2(nodegeo_indices_to_keep, :);
Cinq = C2 - problem_data.r^2;

c_mat = cell(nL, 1);
for j = 1:nL
    c_mat{j} = ak_tqG{j};
end

end


function [dx_normalized,dz] = my_direction_function(c_val, s_val, w_val)

import casadi.*

c = SX.sym('c');
s = SX.sym('s');
w = SX.sym('w');

vec = [c, s, w];

dx_candidate1 = SX([0, -w, s]);
dx_candidate2 = SX([-w, 0, c]);
dx_candidate3 = SX([-s, c, 0]);

abs_vals = [abs(c), abs(s), abs(w)];
minC_val = min(abs_vals);

eps_tol = 1e-10;
range_abs = max(abs_vals) - min(abs_vals);

cond1 = (minC_val <= abs(c) + eps_tol) & (minC_val >= abs(c) - eps_tol);
cond_range_close_to_zero = (range_abs < eps_tol);
final_cond1 = cond1 | cond_range_close_to_zero;

cond2 = (minC_val <= abs(s) + eps_tol) & (minC_val >= abs(s) - eps_tol) & ~final_cond1;
cond3 = (minC_val <= abs(w) + eps_tol) & (minC_val >= abs(w) - eps_tol) & ~final_cond1 & ~cond2;

dx_base = if_else(final_cond1, dx_candidate1, ...
    if_else(cond2, dx_candidate2, if_else(cond3, dx_candidate3, SX.zeros(1,3))));

dx_cross = cross(vec, dx_base);

norm_dx = sqrt(dot(dx_cross, dx_cross));
dx_normalized = dx_cross / norm_dx;

dz = cross(dx_normalized, vec);

F = Function('F', {c, s, w}, {dx_normalized, dz, final_cond1, cond2, cond3});

[dx_val, dz_val, ~, ~, ~] = F(c_val, s_val, w_val);

dx_normalized = full(dx_val);
dz = full(dz_val);

end


function Nid = findNode(node, pos, tol)

if nargin < 3
    tol = 1e-6;
end

Nid = zeros(size(pos,1),1);
dim = size(pos,2);

for i = 1:size(pos,1)
    for j = 1:size(node,1)
        if norm(pos(i,1:dim) - node(j,1:dim)) < tol
            Nid(i) = j;
        end
    end
end

end


function Nd = getNewNd(problem_data, x)

dim = problem_data.dim;
dofgeo = problem_data.dofgeo;
import casadi.*

if isa(problem_data.Nd, 'casadi.SX')
    Nd = problem_data.Nd;
else
    Nd = SX(problem_data.Nd);
end

u_total_size = size(Nd, 1) * dim;
u = SX.zeros(u_total_size, 1);

mc = nnz(problem_data.dofgeo);
u1 = x.u(1:mc);

numeric_indices_for_u = find(dofgeo == 1);
u(numeric_indices_for_u) = u1;

reshaped_u = reshape(u, int64(dim), int64(size(Nd, 1)));
reshaped_u = reshaped_u';

if isa(problem_data.NdFixed, 'casadi.SX')
    NdFixed_SX = problem_data.NdFixed;
else
    NdFixed_SX = SX(problem_data.NdFixed);
end

Nd(:,1:dim) = NdFixed_SX + reshaped_u;

end


function ipmsgeopt_showTruss(x, prob)

p = prob.problem_data;
mc = nnz(p.dofgeo);
n = length(p.Cnlst);
a = x.u(mc+1:mc+n);

dim = p.dim;
Nd_orig = p.Nd;
c = p.Cn;
cl = p.Cnlst;

u_geo_temp = zeros(size(Nd_orig,1) * dim, 1);
u1 = x.u(1:mc);
u_geo_temp(p.dofgeo == 1) = u1;
Nd = p.NdFixed + reshape(u_geo_temp, dim, [])';

Xsize = max(Nd(:,1)) - min(Nd(:,1));
Ysize = max(Nd(:,2)) - min(Nd(:,2));
dimL = [Xsize, Ysize];
if dim == 3
    Zsize = max(Nd(:,3)) - min(Nd(:,3));
    dimL = [dimL, Zsize];
end

maxDeformation = norm(dimL) / 5.0;

u_disp_temp = zeros(size(Nd_orig,1) * dim, 1);
m1 = nnz(p.dof) / p.nL;
u2 = x.u(mc+n+1:mc+n+m1);
u_disp_temp(p.dof(1:length(p.dof)/p.nL) == 1) = u2;

if (max(abs(u_disp_temp)) > 0)
    factor = maxDeformation / max(abs(u_disp_temp));
    Nd(:,1:dim) = Nd(:,1:dim) + factor * 0 * reshape(u_disp_temp, dim, [])';
end

dim0 = dim;
if dim == 2
    Nd(:, dim+1) = 0;
    dim = 3;
end

hold off;
hold on;
axis equal off;

a_norm = a / max(a);
max_a = max(a);
tol = max_a * 1e-3;
rs = dim0 / 12 * max_a;

for i = 1:length(cl)
    cid = cl(i);
    n1 = c(cid, 1);
    n2 = c(cid, 2);
    rt = dim0 / 10 * a_norm(i);
    
    if a(i) > tol
        clr = [0.5, 0.5, 0.5];
        p1 = Nd(n1, 1:dim);
        p2 = Nd(n2, 1:dim);
        ipmsgeopt_drawtube(p1, p2, rt, clr);
        
        [x_s, y_s, z_s] = sphere(50);
        r_s = 0.5 * rs;
        x_s = x_s * r_s; y_s = y_s * r_s; z_s = z_s * r_s;
        surf(x_s + p1(1), y_s + p1(2), z_s + p1(3), 'FaceColor', 'k', 'EdgeColor', 'none');
        surf(x_s + p2(1), y_s + p2(2), z_s + p2(3), 'FaceColor', 'k', 'EdgeColor', 'none');
    end
end

view(dim0);
box on;
axis equal;
axis off;
lighting phong;
material shiny;
camlight('headlight');

end


function ipmsgeopt_drawtube(p1, p2, r, clr)

u = p2 - p1;
t = r / 2 * null(u)';
v = t(1,:);
w = t(2,:);

m = 2;
n = 100;
[s, t_mesh] = meshgrid(linspace(0, 1, m), linspace(0, 2*pi, n));
s = s(:);
t_mesh = t_mesh(:);

p_tube = repmat(p1, m*n, 1) + s * u + cos(t_mesh) * v + sin(t_mesh) * w;
X = reshape(p_tube(:,1), n, m);
Y = reshape(p_tube(:,2), n, m);
Z = reshape(p_tube(:,3), n, m);

hold on
h = surf(X, Y, Z, 'FaceColor', clr, 'FaceAlpha', 1, 'LineStyle', 'none');
h.AmbientStrength = 0.1;
h.DiffuseStrength = 0.1;
h.SpecularStrength = 0.1;
h.SpecularExponent = 1;

hold on

end
