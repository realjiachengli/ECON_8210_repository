% Author: Jiacheng Li
% Date: Oct. 19, 2024
%
% This is the main script for HW 1 of the course ECON 8210. The Miranda-Fackler compecon toolbox is used
% to solve some of the problems.


% ---------------------------------------------------
% 0. Housekeeping
% ---------------------------------------------------

clear all; close all;

addpath('../compecon/CEtools');



% ---------------------------------------------------
% 1. Integration
% ---------------------------------------------------

% parameters
T = 100;
rrho = 0.04; 
llambda = 0.02;
u = @(x) - exp(-x);
a = 0;  % lower limit
b = T;  % upper limit
N = 5000;  % number of nodes

% integrand
fn = @(t) exp(-rrho * t) * u(1 - exp(-llambda * t));


% ----------- midpoint rule -----------
I_mp = midpoint_rule(u, a, b, N);


% ----------- trapezoid rule ----------- 
[x,w] = qnwtrap(N, a, b);   % get the weights and nodes
I_tp = w' * u(x);


% ----------- simpson's rule ----------- 
[x,w] = qnwsimp(N, a, b);   % get the weights and nodes
I_sp = w' * u(x);


% ----------- Monte Carlo ----------- 
rng(0);  % for reproducibility
x_mc = a + (b - a) * rand(N, 1);  % random samples
I_mc = (b - a) * mean(u(x_mc));


% Create a table to compare the results
methods = {'Midpoint Rule', 'Trapezoid Rule', 'Simpson''s Rule', 'Monte Carlo'};
results = [I_mp, I_tp, I_sp, I_mc];

% Display the table
T = table(methods', results', 'VariableNames', {'Method', 'Result'});
disp(T);



% ---------------------------------------------------
% 2. Optimization: basic problem
% ---------------------------------------------------

% objective funciton and gradient
objfn = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;


% gradient and Hessian
grad_objfn = @(x) [
    -400 * x(1) * (x(2) - x(1)^2) + 2 * x(1) - 2;
    200 * (x(2) - x(1)^2)
];

hess_objfn = @(x) [
    1200 * x(1)^2 - 400 * x(2) + 2, -400 * x(1);
    -400 * x(1), 200
];


% initial guess
x0 = [0.2, 0.2]';

% parameters
aalpha = 0.01;  % step size
tol = 1e-6;  % tolerance
max_iter = 10000;  % maximum number of iterations




% -------------- BFGS --------------
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');      % by default, fminunc quasi-Newton uses BFGS method
[x_opt, f_opt] = fminunc(objfn, x0, options);

% display the results
fprintf('\nOptimal solution: x1 = %.4f, x2 = %.4f\n', x_opt(1), x_opt(2));
fprintf('Optimal value of the objective function: %.4f\n', f_opt);


% -------------- steepest descent --------------

% solve using steepest descent
[x_opt, f_opt, x_path] = grad_descent_methods(objfn, grad_objfn, hess_objfn, x0, aalpha, tol, max_iter, 'steepest_descent');

% display the results
fprintf('\nOptimal solution: x1 = %.4f, x2 = %.4f\n', x_opt(1), x_opt(2));
fprintf('Optimal value of the objective function: %.4f\n', f_opt);


% -------------- Newton-Raphson --------------

% solve using Newton-Raphson
[x_opt, f_opt] = grad_descent_methods(objfn, grad_objfn, hess_objfn, x0, aalpha, tol, max_iter, 'Newton-Raphson');

% display the results
fprintf('\nOptimal solution: x1 = %.4f, x2 = %.4f\n', x_opt(1), x_opt(2));
fprintf('Optimal value of the objective function: %.4f\n', f_opt);


% -------------- conjugate descent --------------

% solve using conjugate descent
[x_opt, f_opt] = grad_descent_methods(objfn, grad_objfn, hess_objfn, x0, aalpha, tol, 10000, 'momentum');

% display the results
fprintf('\nOptimal solution: x1 = %.4f, x2 = %.4f\n', x_opt(1), x_opt(2));
fprintf('Optimal value of the objective function: %.4f\n', f_opt);




% ---------------------------------------------------
% 3. Computing Pareto efficient allocations
% ---------------------------------------------------

% -------------- solve a simple problem --------------
n_s = 3; 
m_s = 3; 
llambda_s = [0.5; 0.25; 0.25];
aalpha_s = [1; 1; 1];
oomega_s = [-.5 * ones(3, 1), -.2 * ones(3, 1), -.5 * ones(3, 1)];   % i - agent, j - good
e_s = [20, 0, 0;        % i - agent, j - good
     10, 20, 10;
     0, 0, 20];       


% solve the social planner problem
solve_SP(n_s, m_s, llambda_s, aalpha_s, oomega_s, e_s);



% -------------- m = n = 10 --------------
n = 10; 
m = 10; 
llambda = [10, 1, 1, 1, 1, 1, 1, 1, 1, 1]';
aalpha = rand(m, 1);
oomega = - rand(m, n);   % i - agent, j - good
e = 10 * rand(n, m);     


% solve the social planner problem
solve_SP(n, m, llambda, aalpha, oomega, e);





% ---------------------------------------------------
% 4. Computing Equilibrium allocations
% ---------------------------------------------------
solve_CE(n_s, m_s, aalpha_s, oomega_s, e_s);



% ---------------------------------------------------
% Useful functions
% ---------------------------------------------------

function I = midpoint_rule(f, a, b, N)
    h = (b - a) / N;    % step size
    x_mid = a + (0.5 + (0:N-1)) * h;    % points of evaluation
    I = h * sum(f(x_mid));
end


function [x_opt, f_opt, x_path] = grad_descent_methods(f, grad_f, hess_f, x0, aalpha, tol, max_iter, type)
    % Minimizes a multi-dimensional function using steepest descent.
    % INPUTS
    %    f          : function to minimize
    %    grad_f     : gradient of the function
    %    x0         : initial guess
    %    aalpha      : step size
    %    tol        : tolerance
    %    max_iter   : maximum number of iterations
    %    type       : type of method to use
    % OUTPUTS
    %    x_opt      : optimal sol
    %    f_opt      : optimal value of the function
    %    x_path     : path of x values

    x = x0;
    x_path = zeros(max_iter, length(x0));  % preallocate path matrix
    x_path(1, :) = x;  % initialize path with the initial guess

    d = grad_f(x0) / norm(grad_f(x0));
    beta = 0.5;

    for iter = 1:max_iter
        grad = grad_f(x);       % evaluate gradient
        hess = hess_f(x);       % evaluate Hessian

        if norm(grad) < tol
            x_path = x_path(1:iter, :);  % trim the unused part of the path
            break;
        end

        if strcmp(type, 'steepest_descent')
            d = - aalpha * grad / norm(grad);
            x = x + d;
            ndisplay = 1000;
        elseif strcmp(type, 'Newton-Raphson')
            d = - hess \ grad;
            x = x + d;
            ndisplay = 1;
        elseif strcmp(type, 'momentum')
            d = beta * d - aalpha * grad / norm(grad);
            x = x + d;
            ndisplay = 1000;
        end

        x_path(iter + 1, :) = x;  % store current x in the path

        % display the progress
        if mod(iter, ndisplay) == 0
            fprintf('\nIteration %d: f(x) = %.4f, direction = [%.4f, %.4f], x1 = %.4f, x2 = %.4f\n', ...
                iter, f(x), d(1), grad(2), x(1), x(2));
        end
    end

    x_opt = x;
    f_opt = f(x_opt);
end


function x_opt_matrix = solve_SP(n, m, llambda, aalpha, oomega, e)
    % Solves the social planner problem.
    % INPUTS
    %    n          : Number of agents
    %    m          : Number of goods
    %    llambda    : Elasticity parameters - vector of length n
    %    aalpha     : Weights on goods - sector of length m
    %    oomega     : n x m matrix of omega_j^i values
    %    e          : n x m matrix of endowments
    % OUTPUTS
    %    x_opt_matrix  : Optimal allocations

    % total endowments for each good
    e_total = sum(e, 1);

    % initial guess
    x0 = e_total / n;
    x0 = repmat(x0, n, 1);
    x0 = x0(:);

    % set up constraints to the problem
    lb = ones(n * m, 1) * 1e-5;
    ub = [];

    % equality constraints/resource constraints: sum_i x_j^i = e_total_j for each good j
    Aeq = zeros(m, n * m);
    for j = 1:m
        for i = 1:n
            idx = (j - 1) * n + i; % index for summing over agents
            Aeq(j, idx) = 1;
        end
    end
    beq = e_total';

    % inequality constraints
    A = [];
    b = [];

    objective = @(x) SP_objective(x, n, m, oomega, llambda, aalpha);
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true);

    % solve
    [x_opt, fval, ~, ~] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, [], options);

    % unpack results
    x_opt_matrix = reshape(x_opt, [n, m]);

    disp('Optimal allocations (x_j^i):');
    disp(x_opt_matrix);

    disp('Maximum value of the objective function:');
    disp(-fval);

end


function [f, grad_f] = SP_objective(x, n, m, oomega, llambda, aalpha)
    % Calculates the objective function and its gradient for the SP.
    %
    % INPUTS
    %    x          : Allocations (vector of length n*m)
    %    n          : Number of agents
    %    m          : Number of goods
    %    oomega      : n x m matrix of omega_j^i values
    %    llambda     : Vector of length n
    %    aalpha      : Vector of length m
    %
    % OUTPUTS
    %    f          : Value of the objective function (scalar)
    %    grad_f     : Gradient of the objective function (vector of length n*m)
    % We can compute the gradient because of the nice structure of the problem.

    
        % reshape x into an n x m matrix
        x_matrix = reshape(x, [n, m]);
    
        % Initialize the objective function value
        f = 0;
        
        % initialize gradient
        grad_f_matrix = zeros(n, m);
    
        % compute the objective and gradient
        for i = 1:n
            llambda_i = llambda(i);
            for j = 1:m
                x_ij = x_matrix(i, j);
                oomega_ij = oomega(i, j);
                aalpha_j = aalpha(j);
                
                % % sanity check
                % if x_ij <= 1e-5
                %     f = - 1000;
                %     grad_f = zeros(n*m, 1);
                %     return;
                % end
                
                % objective
                term = aalpha_j * x_ij^(1 + oomega_ij) / (1 + oomega_ij);
                f = f + llambda_i * term;
                
                % gradient
                grad_term = aalpha_j * (1 + oomega_ij) * x_ij^oomega_ij / (1 + oomega_ij);
                grad_f_matrix(i, j) = llambda_i * grad_term;
            end
        end
    
        % Flatten the gradient matrix to a vector
        grad_f = grad_f_matrix(:);
    
        % minimization problem
        f = -f;
        grad_f = -grad_f;
    end


    function xji_sol = solve_CE(n, m, aalpha, oomega, e)
        % Solves the competitive equilibrium
        % INPUTS
        %    n          : Number of agents
        %    m          : Number of goods
        %    aalpha     : Weights on goods - sector of length m
        %    oomega     : n x m matrix of omega_j^i values
        %    e          : n x m matrix of endowments
        % OUTPUTS
        %    xji_sol    : CE allocations

        p0 = ones(m - 1, 1);
        lambda0 = ones(n, 1);
        x0 = [p0; lambda0];

        % objective function
        fun = @(x) equilibrium_conditions(x, n, m, aalpha, oomega, e);
        options = optimoptions('fsolve', 'Display', 'iter');

        % Solve the system
        [x_sol, fval, exitflag, output] = fsolve(fun, x0, options);

        p_sol = [1; x_sol(1:m-1)];
        lambda_sol = x_sol(m:end);

        % Compute the allocations x_j^i
        xji_sol = ((lambda_sol * p_sol') ./ aalpha').^(1 ./ oomega);

        % report
        fprintf('Equilibrium Prices (p_j):\n');
        for j = 1:m
            fprintf('p_%d = %.4f\n', j, p_sol(j));
        end

        fprintf('\nLagrange Multipliers (llambda^i):\n');
        for i = 1:n
            fprintf('llambda^%d = %.4f\n', i, lambda_sol(i));
        end

        fprintf('\nAllocations (x_j^i):\n');
        disp(xji_sol);
    end



    function F = equilibrium_conditions(x, n, m, aalpha, oomega, e)
        % Computes the residuals of the equilibrium conditions for the competitive equilibrium
        %
        % INPUTS:
        %    x      : Vector of variables [p_2; p_3; ...; p_m; llambda^1; llambda^2; ...; llambda^n]
        %             where p_j are prices (excluding the numeraire p_1 = 1) and llambda^i are Lagrange multipliers.
        %    n      : Number of agents.
        %    m      : Number of goods.
        %    aalpha  : Vector of preference weights for each good (length m).
        %    oomega  : n x m matrix of preference parameters omega_j^i for each agent i and good j.
        %    e      : n x m matrix of endowments e_j^i for each agent i and good j.
        %
        % OUTPUT:
        %    F      : Vector of residuals of the equilibrium equations (length m + n - 1).
        %
        % The function computes the residuals of the following equations:
        %    1. Market clearing conditions for each good.
        %    2. Budget constraints for each agent (excluding one redundant constraint).
        
            % Unpack endogenous x
            p = [1; x(1:m-1)];      % p1 = 1
            llambda = x(m:end);
            
            aalpha = aalpha(:);
            
            % get demand function
            llambda_p = llambda * p'; % n x m matrix
            xji = ((llambda_p ./ aalpha').^(1 ./ oomega));
            
            % market clearing conditions 
            market_clearing = sum(xji, 1)' - sum(e, 1)';
            
            % budget constraints 
            BCs = xji * p - e * p;
            budget_constraints = BCs(2:end);
            
            F = [market_clearing; budget_constraints];
        end