function ipnsdp_check_prob_fields(prob)
% IPNSDP_CHECK_PROB_FIELDS Validate the fields of the problem struct `prob`
%
%   Checks whether the input struct `prob` contains only the allowed fields
%   required for the IPNSDP problem setup. It verifies the top-level fields,
%   as well as the subfields in `prob.x0` and `prob.options` if they exist.
%
%   If any unexpected or unknown fields are found, the function throws an
%   error with a message listing the invalid fields.
%
% Input
%   prob : struct
%       The problem definition struct to be checked.

% Throws error if:
%   - `prob` contains any fields not in the allowed list.
%   - `prob.x0` contains unknown fields.
%   - `prob.options` contains unknown fields.
%

% List of allowed top-level fields in prob
allowed_fields = {
    'name', 'problem_data', 'nX', 'dimX', 'nu', 'lbX', 'ubX', 'lbu', 'ubu','nPSDcon', ...
    'x0', 'obj', 'A11', 'A12', 'b1', 'A21', 'A22', 'b2', 'nlcon', ...
    'options', 'solve'
};

% Check for unexpected fields
prob_fields = fieldnames(prob);
extra_fields = setdiff(prob_fields, allowed_fields);

if ~isempty(extra_fields)
    error(['Error: prob contains unexpected fields: ', strjoin(extra_fields, ', ')]);
end

% Check subfields for prob.x0 if it exists
if isfield(prob, 'x0')
    allowed_x0_fields = {'X', 'u'};
    x0_fields = fieldnames(prob.x0);
    extra_x0_fields = setdiff(x0_fields, allowed_x0_fields);
    if ~isempty(extra_x0_fields)
        error(['Error: prob.x0 contains unexpected fields: ', strjoin(extra_x0_fields, ', ')]);
    end
end

% Check subfields for prob.options if it exists
if isfield(prob, 'options')
    allowed_options_fields = {
        'scale', 'opt_tol', 'feas_tol', 'max_iter', 'min_mu', 'sigma', ...
        'tau', 'gamma', 'prob', 'linear_eq_solver', 'stepMode', ...
        'chord_decomp', 'clique_merge', 'density', 'min_block_size', ...
        'rel_block_size', 'linconCheck', 'nonlinconCheck', 'penlpMode', ...
        'penlp', 'penldeg'
    };
    options_fields = fieldnames(prob.options);
    extra_options_fields = setdiff(options_fields, allowed_options_fields);
    if ~isempty(extra_options_fields)
        error(['Error: prob.options contains unexpected fields: ', strjoin(extra_options_fields, ', ')]);
    end
end
end
