% Author: Jiacheng Li
% Date: Oct. 20, 2024
%
% This is the main script for HW 1 of the course ECON 8210. This script contains answer to the value function iteration problem.


% ---------------------------------------------------
% 0. Housekeeping
% ---------------------------------------------------
close all; clear all;
delete(gcp('nocreate'));
% parpool(8);

addpath('..'); 



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

[~, V_ss, y_ss, g_ss] = SSeq(x_ss, params);

disp('Steady state values:');
disp(['Capital (k): ', num2str(k_ss)]);
disp(['Labor (l): ', num2str(l_ss)]);
disp(['Consumption (c): ', num2str(c_ss)]);
disp(['Investment (inv): ', num2str(inv_ss)]);
disp(['Output (y): ', num2str(y_ss)]);
disp(['Gov consumption (g): ', num2str(g_ss)]);
disp(['Value: ', num2str(V_ss)]);




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
% ttau_grid = [0.2, 0.25, 0.3];
% zz_grid = [-0.0673, -0.0336, 0, 0.0336, 0.0673];

ttau_grid = [0.25, 0.25, 0.25];
zz_grid = [-0.0, -0.0, 0, 0.0, 0.0];


exst_cell = {ttau_grid, zz_grid};   % exogenous state grids
prob_cell = {transmat_ttau, transmat_zz};   % transition probabilities
[exstmat, transmat, indlist] = grid.Helpers.makeTotalTransition(exst_cell, prob_cell);

% exogenous envrionemnt
exogenv = struct();
exogenv.exnames = {'ttau', 'zz'};
exogenv.exnpt = size(exstmat, 1);
exogenv.exstmat = exstmat;
exogenv.indlist = indlist;
exogenv.transmat = transmat;
exogenv.meansts = [2, 3];   % neutral states
exogenv.meanstid = 8;   % neutral state index

% state variable grids
kgrid = linspace(0.1, 30, 25);
invgrid = linspace(0.1, 10, 25);
ggrid = grid.TensorGrid({1:exogenv.exnpt, kgrid, invgrid});

% everything into obj
obj = struct('params', params, 'exogenv', exogenv, 'ggrid', ggrid);

% initial guess
Vmat0 = V_ss * ones(obj.ggrid.Npt, 1);
Pmat0 = [c_ss, l_ss] .* ones(obj.ggrid.Npt, 1);

maxIter = 1000;
tol = 1e-5;
Verr = 100;

Vmat = Vmat0;
Pmat = Pmat0;

iter = 0;



% value function iteration
while iter < maxIter && Verr > tol

    if mod(iter, 10) == 0
        accel = 0;
    else
        accel = 1;
    end
    
    [Vmat, Pmat] = updateV(Vmat, Pmat, obj, accel);

    Verr = max(abs(Vmat - Vmat0));
    disp(['Iteration ', num2str(iter), ', Vf error: ', num2str(Verr), ', Accelerate: ', num2str(accel)]);
    
    Vmat0 = Vmat;
    iter = iter + 1;
end

% plot the value function
Vf = grid.LinearInterpFunction(obj.ggrid, Vmat);
Pf = grid.LinearInterpFunction(obj.ggrid, Pmat);


grid_plot = [exogenv.meanstid * ones(length(kgrid), 1), kgrid', inv_ss * ones(length(kgrid), 1)];

% Plot the value function
figure;
plot(kgrid, Vf.evaluateAt(grid_plot), 'LineWidth', 2);
title('Value Function over Capital');
xlabel('Capital (k)');
ylabel('Value Function (V)');
grid on;

% Plot the policy functions
pol_plot = Pf.evaluateAt(grid_plot);

figure;
subplot(2, 1, 1);
plot(kgrid, pol_plot(1, :), 'LineWidth', 2);
title('Policy Function for Consumption over Capital');
xlabel('Capital (k)');
ylabel('Consumption (c)');
grid on;

subplot(2, 1, 2);
plot(kgrid, pol_plot(2, :), 'LineWidth', 2);
title('Policy Function for Labor over Capital');
xlabel('Capital (k)');
ylabel('Labor (l)');
grid on;







% ---------------------------------------------------
% Useful functions
% ---------------------------------------------------

function [res, V, y, g] = SSeq(x, params)
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
    eeta   = params.eeta;      % utility weight on government consumption

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
    res(2) = l - (1 / c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * (1 - aalpha) * y;
    res(3) = c + inv - ((1 - ttauSS) * (1 - aalpha) + aalpha) * y;
    res(4) = (1 / c) - bbeta * (1/c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * aalpha * y;

    % steady state values
    g = ttauSS * (1 - aalpha) * y;
    V = log(c) + eeta * log(g) - l^2 / 2;
end


function [Vmat, Pmat] = updateV(Vmat, Pmat, obj, accel)
    % Each iteration in fixed grid VFI
    % INPUTs
    %     Vf: interpolated value function
    %     ggrid: a structure with grids
    %     exogenv: a structure containing the exogenous environment
    %     params: a structure containing the model parameters
    % OUTPUTs
    %     Vf: updated value function

    % Unpack the parameters
    aalpha = obj.params.aalpha;
    bbeta = obj.params.bbeta;
    eeta = obj.params.eeta;

    % unpack grid points
    gridmat = obj.ggrid.Pointmat;
    Npt = obj.ggrid.Npt;

    % interpolate the value function
    Vf = grid.LinearInterpFunction(obj.ggrid, Vmat);

    for ii = 1:Npt

        if ~accel
            % unpack the state
            state = gridmat(ii, :);
            ttau = obj.exogenv.exstmat(gridmat(ii, 1), 1);
            zz = obj.exogenv.exstmat(gridmat(ii, 1), 2);
    
            k = state(2);
    
            % upper bound on consumption
            cmax = ((1 - ttau) * (1 - aalpha) + aalpha) * exp(zz) * k^aalpha / (1 - aalpha)^((aalpha - 1) / 2);
            obj2min = @(c) - pVc(c, gridmat(ii, :), Vf, obj);
    
            % minimize - value
            opts = optimset('Display', 'off');
            Pmat(ii, 1) = fminbnd(obj2min, 1e-6, cmax, opts);
        end

        [pV, Pmat(ii, 2), g] = pVc(Pmat(ii, 1), gridmat(ii, :), Vf, obj);

        Vmat(ii) = pV + (1 - bbeta) * eeta * log(g);
    end

end


function [val, lopt, g] = pVc(c, state, Vf, obj)

    % unpack params and states
    aalpha = obj.params.aalpha;
    ddelta = obj.params.ddelta;
    pphi = obj.params.pphi;
    eeta = obj.params.eeta;
    bbeta = obj.params.bbeta;

    ttau = obj.exogenv.exstmat(state(1), 1);
    zz = obj.exogenv.exstmat(state(1), 2);

    k = state(2);
    inv = state(3);

    % optimal labor
    lopt = ((1 - aalpha) * ((1 - ttau) * (1 - aalpha) + aalpha) * exp(zz) * k^aalpha / c)^(1 / (1 + aalpha));
    y = exp(zz) * k^aalpha * lopt^(1 - aalpha);
    income = ((1 - ttau) * (1 - aalpha) + aalpha) * y;
    g = ttau * (1 - aalpha) * y;

    % next state
    inv_next = income - c;
    k_next = (1 - ddelta) * k + (1 - pphi * (inv_next / inv - 1)^2) * inv_next;

    state_next = [(1:obj.exogenv.exnpt)', [k_next, inv_next] .* ones(obj.exogenv.exnpt, 1)];
    Vf_next = Vf.evaluateAt(state_next);
    expVf_next = obj.exogenv.transmat(state(1), :) * Vf_next';

    % present value
    val = (1 - bbeta) * (log(c) - lopt^2 / 2) + bbeta * expVf_next;
end