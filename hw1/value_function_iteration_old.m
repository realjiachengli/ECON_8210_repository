% Author: Jiacheng Li
% Date: Oct. 20, 2024
%
% This is the main script for HW 1 of the course ECON 8210. The Miranda-Fackler compecon toolbox is used
% to solve some of the problems.



% ---------------------------------------------------
% 0. Housekeeping
% ---------------------------------------------------
close all; clear all;

addpath('..');  % to make use of +grid and +simul. The utilities are developed in Elenev et al (2022)



% ---------------------------------------------------
% 1. Compute the steady state
% ---------------------------------------------------

% Define the parameters
params.aalpha = 0.33;
params.bbeta = 0.97;
params.ddelta = 0.1;
params.pphi = 0.05;
params.eeta = 0.2;
params.ttauSS = 0.25;
params.zzSS = 0;


% Initial guess for the endogenous variables [k, l, c, inv]
x0 = [1, 0.5, 0.5, 0.5]';

% Solve the steady state equations using fsolve
options = optimoptions('fsolve', 'Display', 'iter');
[x_ss, fval, exitflag] = fsolve(@(x) SSeq(x, params), x0, options);

% Display the results
k_ss = x_ss(1);
l_ss = x_ss(2);
c_ss = x_ss(3);
inv_ss = x_ss(4);

disp('Steady state values:');
disp(['Capital (k): ', num2str(k_ss)]);
disp(['Labor (l): ', num2str(l_ss)]);
disp(['Consumption (c): ', num2str(c_ss)]);
disp(['Investment (inv): ', num2str(inv_ss)]);




% ---------------------------------------------------
% 2. VFI with fixed grid
% ---------------------------------------------------

% shocks
transmat_ttau = [0.9, 0.1, 0; 0.05, 0.9, 0.05; 0, 0.1, 0.9];
transmat_zz = [
        0.9727 0.0273 0      0      0;
        0.0041 0.9806 0.0153 0      0;
        0      0.0082 0.9836 0.0082 0;
        0      0      0.0153 0.9806 0.0041;
        0      0      0      0.0273 0.9727
    ];
ttau_grid = [0.2, 0.25, 0.3];
zz_grid = [-0.0673, -0.0336, 0, 0.0336, 0.0673];

exst_cell = {ttau_grid, zz_grid};   % exogenous state grids
prob_cell = {transmat_ttau, transmat_zz};   % transition probabilities
[exstmat, transmat, indlist] = grid.Helpers.makeTotalTransition(exst_cell, prob_cell);

% exogenous envrionemnt
exogenv = struct;
exogenv.exnames = {'ttau', 'zz'};
exogenv.exnpts = size(exstmat, 1);
exogenv.exstmat = exstmat;
exogenv.transmat = transmat;

% endogenous state grid
k_grid = linspace(0.7 * k_ss, 1.3 * k_ss, 20);
inv_grid = linspace(inv_ss - 0.5 * k_ss, inv_ss + 0.5 * k_ss, 10);

enst_cell = {k_grid, inv_grid};   % endogenous state grids
st_cell = [{1:exogenv.exnpts}, enst_cell];   % total state grids

mygrid = grid.TensorGrid(st_cell);

% put everything useful together in an object
obj = struct;
obj.params = params;
obj.exogenv = exogenv;
obj.grid = mygrid;


solve_VFI_fixed_grid(x_ss, obj);



% ---------------------------------------------------
% 3. VFI with endogenous grid
% ---------------------------------------------------






% ---------------------------------------------------
% 4. Comparison of grids
% ---------------------------------------------------






% ---------------------------------------------------
% 5. Switching between PFI and VFI
% ---------------------------------------------------






% ---------------------------------------------------
% 6. Multigrid
% ---------------------------------------------------






% ---------------------------------------------------
% 7. Stochastic grid
% ---------------------------------------------------




% ---------------------------------------------------
% Useful functions
% ---------------------------------------------------

function res = SSeq(x, params)
    % Steady state characterization equations
    % INPUTs
    %     x: a vector of endogenous variables
    %     params: a structure containing the model parameters
    % OUTPUTs
    %     res: a vector of residuals

    % Unpack the parameters
    aalpha  = params.aalpha;    % capital share
    bbeta   = params.bbeta;     % discount factor
    ddelta  = params.ddelta;    % depreciation rate
    pphi    = params.pphi;    % capital adjustment cost
    eeta    = params.eeta;    % gov consumption weight

    ttauSS  = params.ttauSS;    % steady state tax on labor
    zzSS    = params.zzSS;      % steady state productivity
    

    % Unpack the endogenous variables
    k = x(1);   % capital
    l = x(2);   % labor
    c = x(3);   % consumption
    inv = x(4);  % investment

    
    % intermediate variables
    y = exp(zzSS) * k^aalpha * l^(1-aalpha);    % output
    
    % Compute the steady state residuals
    res = ones(4, 1);

    res(1) = ddelta * k - inv;
    res(2) = eeta * (1 - aalpha) * (1 / l) - (l - (1 / c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * (1 - aalpha) * y);
    res(3) = c + inv - ((1 - ttauSS) * (1 - aalpha) + aalpha) * y;
    res(4) = (1 / c) - bbeta * (eeta * aalpha * (1 / k) + (1/c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * aalpha * y);


end


function transpath = simul_IRF(obj)
end


function [Vf, Pf] = solve_VFI_fixed_grid(sol_ss, obj)
    % Solve the model using value function iteration with fixed grid
    % INPUTs
    %     params: a structure containing the model parameters
    % OUTPUTs
    %     Vf: a 1D array of value functions [V]
    %     Pf: a 2D array of policy functions [c, l]

    pointmat = obj.grid.Pointmat;
    params = obj.params;
    npts = size(pointmat, 1);

    % unpack the steady state and parameters
    k_ss = sol_ss(1);
    l_ss = sol_ss(2);
    c_ss = sol_ss(3);
    inv_ss = sol_ss(4);

    ttauSS = params.ttauSS;
    zzSS = params.zzSS;
    aalpha = params.aalpha;
    eeta = params.eeta;

    % initialize the value function
    Vmatguess = ones(npts, 1) * ...
        (log(c_ss) + eeta * log(ttauSS * (1 - aalpha) * exp(zzSS) * k_ss^aalpha * l_ss^(1 - aalpha)) - (l_ss^2) / 2);

    Vf = grid.LinearInterpFunction(obj.grid, Vmatguess);

    x0 = [c_ss, l_ss] .* ones(npts, 1);
    
    % set hyperparameters
    maxiter = 1000;
    Vf_tol = 1e-6;

    Vmat0 = Vmatguess;

    % start counter
    iter = 0;

    % value function iteration
    while iter < maxiter
        iter = iter + 1;

        % define the present value to maximize
        objfunc = @(sol) - pv(sol, obj, Vf);

        % solve for the optimal policy
        options = optimoptions('fminunc', 'Display', 'iter');
        x_sol = fminunc(objfunc, x0, options);

        Vmat = - objfunc(x_sol);
        dist = max(abs(Vmat0 - Vmatguess));
        if dist < Vf_tol
            disp('Vmat converged.')
            break;
        end

        disp(['Iteration: ', num2str(iter), ' max dist: ', num2str(dist)]);

        % update the value function
        Vf = Vf.fitTo(Vmat);

    end

end


function V = pv(sol, obj, Vf)
    % Compute the present value of a given sol
    % INPUTs
    %     sol: a 2D array of endogenous variables [c, l]
    %     obj: a structure containing the model parameters, grids and exogenous environment
    % OUTPUTs
    %      V: a 1D array of present values

    % unpack exogenous environment
    exstmat = obj.exogenv.exstmat;
    transmat = obj.exogenv.transmat;
    pointmat = obj.grid.Pointmat;
    npts = size(pointmat, 1);
    exnpts = size(exstmat, 1);

    % unpack grid points
    exstid = pointmat(:, 1);
    k = pointmat(:, 2);
    inv = pointmat(:, 3);
    ttau = exstmat(exstid, 1);
    zz = exstmat(exstid, 2);

    % unpack variables to solve for
    c = sol(:, 1);
    l = sol(:, 2);

    % unpack parameters
    aalpha = obj.params.aalpha;
    bbeta = obj.params.bbeta;
    ddelta = obj.params.ddelta;
    pphi = obj.params.pphi;
    eeta = obj.params.eeta;

    % current period value
    curr_V = (1 - bbeta) * (log(c) + eeta * log(ttau .* (1 - aalpha) .* exp(zz) .* k.^aalpha .* l.^(1 - aalpha)) - (l.^2) / 2);

    % state transitions
    inv_next = ((1 - ttau) * (1 - aalpha) + aalpha) .* exp(zz) .* k.^aalpha .* l.^(1 - aalpha) - c;
    k_next = (1 - ddelta) * k + (1 - pphi * ((inv_next ./ inv) - 1).^2) .* inv_next;

    % compute expected continuation value
    nextpr = transmat(exstid, :);

    E = @(V) sum(nextpr .* V, 2);

    nextst = [repmat((1:exnpts), npts, 1), kron([inv_next, k_next], ones(1, exnpts))];   % [exstid, k_next, inv_next] all npts by exnpts
    nextstmat = reshape(nextst, npts, exnpts, 3);   % reshape to npts by exnpts by nstates

    % reshape again to do evaluation
    nextst = reshape(nextstmat, npts * exnpts, 3);
    V_next = Vf.evaluateAt(nextst);
    V_next = reshape(V_next, npts, exnpts);  % reshape back to npts by exnpts

    % compute expected value
    expV = E(V_next);

    % total present value
    V = curr_V + bbeta * expV;

end
