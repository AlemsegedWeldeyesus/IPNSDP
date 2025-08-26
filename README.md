# IPNSDP: A Solver for Nonlinear Semidefinite Programming with Automatic Differentiation and Chordal Decomposition

IPNSDP is a MATLAB-based solver for **Nonlinear Semidefinite Programming (NSDP)** problems. It uses an interior-point method to handle symmetric matrix variables (`X_k`), scalar variables (`u`), and both linear and nonlinear constraints, including positive semidefinite (PSD) constraints. To improve scalability for large PSD constraints, IPNSDP supports **chordal decomposition**, which exploits structured sparsity by decomposing large semidefinite constraints into smaller ones. The solver integrates with the **CasADi** framework for efficient computation of derivatives.

---

## Generic NSDP Formulation (Plain-text version)

The general form of an NSDP problem solved by IPNSDP is:

```
minimize f(X, u)
subject to
  sum_k <A11_i^(k), X_k> + A12_i * u = b1_i,  i=1..n1
  sum_k <A21_j^(k), X_k> + A22_j * u <= b2_j, j=1..n2
  c1(X,u) = 0
  c2(X,u) <= 0
  A_s(X,u) >= 0, s=1..n_psd
  u_lower <= u <= u_upper
  X_lower_k * I <= X_k <= X_upper_k * I, k=1..mX
```

---

## Features

- **Problem Types**: Supports symmetric matrix variables, scalar variables, and mixed linear, nonlinear, and PSD constraints  
- **CasADi Integration**: Automatic differentiation with CasADi for efficient derivative computation  
- **Customization**: Flexible options via `prob.options` for tolerances, step strategies, and decomposition methods  
- **Chordal Decomposition**: Exploits structured sparsity in PSD constraints with chordal decomposition and merging  

---

## Installation

1. **CasADi**: Download and install from [CasADi downloads](https://web.casadi.org/get/)  
2. **IPNSDP**: Download from [GitHub](https://github.com/AlemsegedWeldeyesus/IPNSDP)  
3. **Setup**: In MATLAB, navigate to the IPNSDP directory and run:

```matlab
install()
```

This adds the solver to your MATLAB path including all subfolders.

4. **Verify**: Run the test script:

```matlab
ipnsdp_test()
```

A successful run confirms IPNSDP is installed correctly and ready to use.

---

## Problem Setup

IPNSDP uses a structured MATLAB template, `ipnsdp_input_template.m`, to define problems via a structured `prob` object.

### Minimal Required Fields

- `prob.name`: Problem name (string, optional)  
- `prob.problem_data`: User-defined struct for passing data (optional)  

#### Matrix Variables

- `prob.nX`: Number of matrix variables `X1, ..., XmX`. Set to 0 if none  
- `prob.dimX`: Vector `[d_1, ..., d_mX]` of matrix dimensions. `[]` if none  
- `prob.lbX`: Vector `[rho_1, ..., rho_mX]` for lower bounds `Xk >= rho_k * I`. Use `-Inf` for selectively unbounded entries. Defaults used if not provided  
- `prob.ubX`: Vector `[rho_bar1, ..., rho_barmX]` for upper bounds `Xk <= rho_bar_k * I`. Use `Inf` for selectively unbounded entries. Defaults used if not provided  

#### Scalar Variables

- `prob.nu`: Number of scalar variables `u in R^n`. 0 if none  
- `prob.lbu`: Lower bounds `u >= u_lower`. `-Inf` for selectively unbounded entries. Defaults used if not provided  
- `prob.ubu`: Upper bounds `u <= u_upper`. `Inf` for selectively unbounded entries. Defaults used if not provided  

#### PSD Constraints

- `prob.nPSDcon`: Number of nonlinear PSD constraints `A_s(X,u) >= 0`. 0 if none  

#### Linear Constraints

- `prob.A11`: Coefficient matrix for `vec(X)` in linear equality constraints. `[]` if none  
- `prob.A12`: Coefficient matrix for `u` in linear equality constraints. `[]` if none  
- `prob.A21`: Coefficient matrix for `vec(X)` in linear inequality constraints. `[]` if none  
- `prob.A22`: Coefficient matrix for `u` in linear inequality constraints. `[]` if none  

#### Function Handles

- `prob.obj`: Objective function handle, e.g., `@(x) objective(x, prob)`  
- `prob.nlcon`: Nonlinear constraints function handle, returning `[eq, ineq, psd]`  
  - `eq`: Nonlinear equality vector. `[]` if none  
  - `ineq`: Nonlinear inequality vector. `[]` if none  
  - `psd`: Cell array of symmetric matrices for PSD constraints. `{}` if none  

#### Initial Point (Optional)

- `prob.x0.X`: Cell array of initial matrix points, e.g., `x0.X{1} = eye(d1)`. Defaults if not provided  
- `prob.x0.u`: Vector of initial scalar points. Defaults if not provided  

---

## Solver Options

Some options in `prob.options`:

- `chord_decomp`: `'yes'` / `'no'`  
- `clique_merge`: `'yes'` / `'no'`  
- `opt_tol`, `feas_tol`: Stopping tolerances  
- `max_iter`: Maximum iterations  
- `stepMode`: `'joint'` or `'separate'`  
- `linconCheckN`: `'yes'` / `'no'` for linear dependency check  (Numeric)
- `nonlinconCheckN`: `'yes'` / `'no'` for nonlinear dependency check (Numeric)
---

## Output

The solver returns two main outputs: `x` and `info`.

- `x`: Final solution including `x.X`, `x.u`, duals, and slacks  
- `info`: Diagnostics including final objective, number of iterations, CPU time, and exit message  

---

## Documentation

- For clarity of the mathematical formulation, see **ipnsdp.pdf**.  
- For further details on problem setup and solver options, refer to **ipnsdp_input_template.m**.  
- Several illustrative examples are also included in the **examples** directory.  

## License

See the license file.

## Authors

Alemseged Weldeyesus and Miguel F. Anjos

**Correspondence:** a.weldeyesus@ed.ac.uk  
