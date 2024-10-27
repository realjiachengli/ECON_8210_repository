% Author: Jiacheng Li
% Date: Oct. 20, 2024
%
% This is the main script for HW 1 of the course ECON 8210. The Miranda-Fackler compecon toolbox is used
% to solve some of the problems.



% ---------------------------------------------------
% 0. Housekeeping
% ---------------------------------------------------
close all; clear all;
delete(gcp('nocreate'));

addpath('..');  % to make use of +grid and +simul. The utilities are developed in Elenev et al (2022)



% ---------------------------------------------------
% 1. Compute the steady state
% ---------------------------------------------------

% Define the parameters
params.aalpha = 0.33;
params.bbeta = 0.97;
params.ddelta = 0.1;
params.pphi = 0.0;
params.eeta = 0.0;
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
y_ss = exp(params.zzSS) * k_ss^params.aalpha * l_ss^(1 - params.aalpha);
g_ss = params.ttauSS * (1 - params.aalpha) * y_ss;

disp('Steady state values:');
disp(['Capital (k): ', num2str(k_ss)]);
disp(['Labor (l): ', num2str(l_ss)]);
disp(['Consumption (c): ', num2str(c_ss)]);
disp(['Investment (inv): ', num2str(inv_ss)]);
disp(['Output (y): ', num2str(y_ss)]);
disp(['Gov consumption (g): ', num2str(g_ss)]);





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
exogenv.indlist = indlist;
exogenv.transmat = transmat;

% endogenous state grid
k_grid = linspace(0.1,20, 25);
inv_grid = linspace(0.1,10, 25);
% k_grid = linspace(0.7 * k_ss, 1.3 * k_ss, 25);
% inv_grid = linspace(inv_ss - 0.5 * k_ss, inv_ss + 0.5 * k_ss, 5);

enst_cell = {k_grid, inv_grid};   % endogenous state grids
st_cell = [{1:exogenv.exnpts}, enst_cell];   % total state grids

mygrid = grid.TensorGrid(st_cell);

% put everything useful together in an object
obj = struct;
obj.params = params;
obj.exogenv = exogenv;
obj.grid = mygrid;
obj.SSvals = [k_ss, l_ss, c_ss, inv_ss, y_ss, g_ss];
obj.accelerate = 1;     % turn on acceleration steps (keep policy fixed for 10 iterations)

% delete(gcp('nocreate'))
% parpool('local', 6);

tic; 
[Vf_solved, Pf_solved] = solve_VFI_fixed_grid(obj);
toc;

obj.Vf = Vf_solved;
obj.Pf = Pf_solved;

% plot the value and policy function along the dimension k
kgrid_plot = linspace(0.7 * k_ss, 1.3 * k_ss, 100)';
Plotgridmat = [8 * ones(100, 1), kgrid_plot, inv_ss * ones(100, 1)];
Vf_grid = obj.Vf.evaluateAt(Plotgridmat);
Pf_grid = obj.Pf.evaluateAt(Plotgridmat);

% Plot the value function
figure;
plot(kgrid_plot, Vf_grid);
title('Value Function');
xlabel('Capital (k)');
ylabel('Value');

% Plot the policy function for consumption
figure;
plot(kgrid_plot, Pf_grid);
title('Policy Function for Consumption');
xlabel('Capital (k)');
ylabel('Consumption (c)');

% Plot the policy function for labor
figure;
plot(kgrid_plot, Pf_grid);
title('Policy Function for Labor');
xlabel('Capital (k)');
ylabel('Labor (l)');

% plot IRF
NT = 100;
shock_id = [2, 1];

simul_IRF(shock_id, NT, obj);


% ---------------------------------------------------
% 3. VFI with endogenous grid
% ---------------------------------------------------
solve_VFI_endogrid(obj);





% ---------------------------------------------------
% 4. Comparison of grids
% ---------------------------------------------------






% ---------------------------------------------------
% 5. Switching between PFI and VFI
% ---------------------------------------------------

% Now turn it off and let's compare speed
obj.accelerate = 0; 

tic; 
[Vf_solved, Pf_solved] = solve_VFI_fixed_grid(obj);
toc;





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
    res(2) = l - (1 / c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * (1 - aalpha) * y;
    res(3) = c + inv - ((1 - ttauSS) * (1 - aalpha) + aalpha) * y;
    res(4) = (1 / c) - bbeta * (1/c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * aalpha * y;


end


function transpath = simul_IRF(shock_id, NT, obj)
    % Simulate the impulse response function
    % INPUTs
    %     sol_ss: a 1D array of steady state values [k, l, c, inv]
    %     shock_id: a 1D array of shock id [ttau_id, zz_id]
    %     NT: a scalar of number of periods to simulate
    %     obj: a structure containing the model parameters, grids and exogenous environment
    % OUTPUTs
    %    transpath: a 2D array of simulated path [exstid, k, inv, c, l, V, ttau, zz]

    % Unpack the steady state values
    k_ss = obj.SSvals(1);
    l_ss = obj.SSvals(2);
    c_ss = obj.SSvals(3);
    inv_ss = obj.SSvals(4);


    % ttauSS_id = 2; 
    % zzSS_id = 3;

    % get shock
    shock_st = find(obj.exogenv.indlist(:, 1) == shock_id(1) & obj.exogenv.indlist(:, 2) == shock_id(2));
    shockSS_id = find(obj.exogenv.indlist(:, 1) == 2 & obj.exogenv.indlist(:, 2) == 3);

    % simulate the path
    transpath = zeros(NT, 7);   % [exstid, k, inv, c, l, V, ttau, zz]
    transpath(:, 1) = shockSS_id;
    transpath(:, 2) = k_ss;
    transpath(:, 3) = inv_ss;
    transpath(:, 4) = c_ss;
    transpath(:, 5) = l_ss;
    transpath(:, 6) = 0;    % value function

    transpath(:, 7) = obj.params.ttauSS;
    transpath(:, 8) = obj.params.zzSS;


    for tt = 1:NT
        if tt == 1  % initial state
            transpath(tt, 1) = shock_st;
            transpath(tt, 7) = obj.exogenv.exstmat(shock_st, 1);
            transpath(tt, 8) = obj.exogenv.exstmat(shock_st, 2);
        end

        curr_st = transpath(tt, 1:3);

        % policy path
        pol_sol = obj.Pf.evaluateAt(curr_st);
        transpath(tt, 4) = pol_sol(1);
        transpath(tt, 5) = pol_sol(2);
        transpath(tt, 6) = obj.Vf.evaluateAt(curr_st);

        % state transition
        enst_next = enst_transition(curr_st, pol_sol, obj.params, obj.exogenv);
        transpath(tt + 1, 2:3) = enst_next;
    end

    % plot transpath in subplots
    figure;
    subplot(3, 2, 1);
    plot(transpath(:, 2));
    title('Capital (k)');
    xlabel('Time');
    ylabel('k');

    subplot(3, 2, 2);
    plot(transpath(:, 3));
    title('Investment (inv)');
    xlabel('Time');
    ylabel('inv');

    subplot(3, 2, 3);
    plot(transpath(:, 4));
    title('Consumption (c)');
    xlabel('Time');
    ylabel('c');

    subplot(3, 2, 4);
    plot(transpath(:, 5));
    title('Labor (l)');
    xlabel('Time');
    ylabel('l');

    subplot(3, 2, 5);
    plot(transpath(:, 6));
    title('Value Function (V)');
    xlabel('Time');
    ylabel('V');

    subplot(3, 2, 6);
    plot(transpath(:, 7));
    hold on;
    plot(transpath(:, 8));
    title('Shocks (ttau and zz)');
    xlabel('Time');
    ylabel('Shock values');
    legend('ttau', 'zz');
    hold off;



end


function [Vf, Pf] = solve_VFI_fixed_grid(obj)
    % Solve the model using value function iteration with fixed grid
    % INPUTs
    %     sol_ss: a 1D array of steady state values [k, l, c, inv]
    %     obj: a structure containing the model parameters, grids and exogenous environment
    % OUTPUTs
    %     Vf: a 1D array of value functions [V]
    %     Pf: a 2D array of policy functions [c, l]

    pointmat = obj.grid.Pointmat;
    params = obj.params;
    npts = size(pointmat, 1);
    exogenv = obj.exogenv;

    % useful parameters
    aalpha = params.aalpha;
    bbeta = params.bbeta;
    eeta = params.eeta;

    % unpack the steady state and parameters
    k_ss = obj.SSvals(1);
    l_ss = obj.SSvals(2);
    c_ss = obj.SSvals(3);
    inv_ss = obj.SSvals(4);
    g_ss = obj.SSvals(6);

    % initialize the value function
    Vpts0 = ones(npts, 1) * ...
        (log(c_ss) + params.eeta * log(g_ss) - (l_ss^2) / 2);

    Vf = grid.LinearInterpFunction(obj.grid, Vpts0);
    x0 = [0.1, 0.5] .* ones(npts, 1);
    
    % set hyperparameters
    maxiter = 1000;
    Vf_tol = 1e-6;

    % value - policy update
    Vpts = Vpts0;
    xsol = x0;

    % start counter
    iter = 0;

    % value function iteration
    while iter < maxiter

        % parfor ii = 1:npts
        for ii = 1:npts

            curr_st = pointmat(ii, :);
            exst = curr_st(1);
            k = curr_st(2);


            % define the present value to maximize
            objfunc = @(sol) - pv_ng(sol, curr_st, Vf, exogenv, params);

            if mod(iter, 10) == 0       % only solve for optimal policy in these steps
                opts = optimoptions('fmincon', 'Display', 'off');
                xsol(ii, :) = fmincon(objfunc, xsol(ii, :), [], [], [], [], [0.01, 0.01], [inf, inf], [], opts);
            end

            % compute government consumption
            choice = xsol(ii, :);
            l = choice(2);
            
            ttau = exogenv.exstmat(exst, 1);
            zz = exogenv.exstmat(exst, 2);
            g = ttau * (1 - aalpha) * exp(zz) * k^aalpha * l^(1 - aalpha);

            % update the value function vector
            Vpts(ii) = - objfunc(xsol(ii, :)) + (1 - bbeta) * eeta * log(g);
        end

        % check convergence
        dist = max(abs(Vpts - Vpts0));

        if dist < Vf_tol
            disp('Vpts converged.')
            break;
        end

        disp(['Iteration: ', num2str(iter), ' max dist: ', num2str(dist), ' solve optimal: ', num2str(mod(iter, 10) == 0 )]);

        % update the value function
        Vf = Vf.fitTo(Vpts);
        Vpts0 = Vpts;

        iter = iter + 1;
    end

    % get policy function as well
    Pf = grid.LinearInterpFunction(obj.grid, xsol);

end




function V = pv_ng(sol, curr_st, Vf, exogenv, params)
    % Compute the present value of a given sol
    % INPUTs
    %     sol: a 2D array of endogenous variables [c, l]
    %     curr_st: a 1D array of current state [exstid, k, inv]
    %     Vf: a function to evaluate the value function
    %     exogenv: a structure containing the exogenous environment
    %     params: a structure containing the model parameters
    % OUTPUTs
    %      V: a 1D array of present values

    % unpack state
    shockid = curr_st(1);
    k =  curr_st(2);
    inv =  curr_st(3);

    ttau = exogenv.exstmat(shockid, 1);
    zz = exogenv.exstmat(shockid, 2);
    nextpr = exogenv.transmat(shockid, :);

    % unpack parameters
    aalpha = params.aalpha;
    bbeta = params.bbeta;
    ddelta = params.ddelta;
    pphi = params.pphi;

    % unpack policy
    inv_next = sol(1);
    l = sol(2);

    % current period value
    c = ((1 - ttau) * (1 - aalpha) + aalpha) * exp(zz) * k^aalpha * l^(1 - aalpha) - inv_next;
    curr_V = (1 - bbeta) * (log(c) - (l^2) / 2);

    % state transitions
    k_next = (1 - ddelta) * k + inv_next - pphi * ((inv_next / inv) - 1)^2 * inv_next;

    % compute expected continuation value
    nextst = [(1:exogenv.exnpts)', [k_next, inv_next] .* ones(exogenv.exnpts, 1)];   % [exstid, k_next, inv_next]

    % do evaluation
    nextV = Vf.evaluateAt(nextst);

    % compute expected value
    expV = sum(nextpr .* nextV, 2);

    % total present value
    V = curr_V + params.bbeta * expV;
end



function [] = solve_VFI_endogrid(obj)

    pointmat = obj.grid.Pointmat;
    params = obj.params;
    npts = size(pointmat, 1);

    % unpack the steady state and parameters
    k_ss = obj.SSvals(1);
    l_ss = obj.SSvals(2);
    c_ss = obj.SSvals(3);
    inv_ss = obj.SSvals(4);
    g_ss = obj.SSvals(6);

    % initialize the value function
    Vpts = ones(npts, 1) * ...
        (log(c_ss) + params.eeta * log(g_ss) - (l_ss^2) / 2);

    Vf = grid.LinearInterpFunction(obj.grid, Vpts);


    % value function iteration
    while iter < maxiter
        iter = iter + 1;

        [kpts, Vf_new, Vpts_new] = iteration_EGM(obj, Vf, Vpts);

        % check convergence
        dist = max(abs(Vpts_new - Vpts));

        if dist < Vf_tol
            disp('Vpts converged.')
            break;
        end

        disp(['Iteration: ', num2str(iter), ' max dist: ', num2str(dist), ' solve optimal: ', num2str(mod(iter, 10) == 0 )]);

    end


end


function [kpts, Vf_new, Vpts_new] = iteration_EGM(obj, Vf, Vpts)
    % Each iteration of the EGM algorithm
    % INPUTs
    %     obj: a structure containing the model parameters, grids and exogenous environment
    %     Vpts: a 1D array of value functions on the base gridpoints
    % OUTPUTs
    %     kpts: a 1D array of endogenous gridpoints
    %     Vf: the iterpolated value function object
    %     Vpts: a 1D array of the new value functions on the base gridpoints

    % unpack the grid
    pointmat = obj.grid.Pointmat;
    npts = size(pointmat, 1);
    l_ss = obj.SSvals(2);

    for ii = 1:npts

        % goal: evaluate the new value function at each possible i'
        curr_st = pointmat(ii, :);
        curr_exst = curr_st(1);
        curr_k = curr_st(2);
        curr_inv = curr_st(3);

        % shocks
        ttau = obj.exogenv.exstmat(curr_exst, 1);
        zz = obj.exogenv.exstmat(curr_exst, 2);

        % evaluate some current state variables
        curr_Y = ((1 - ttau) * (1 - aalpha) + aalpha) * exp(zz) * curr_k^aalpha * l_ss^(1 - aalpha);
        
        % expected value for each i'
        inv_next = pointmat(:, 3);

        tempV = tempVfunc(inv_next, obj, Vf, curr_st);
        dtempVdinv = (tempV - tempVfunc(inv_next, obj, Vf+1e-5, curr_st)) / 1e-5;  % numerical derivative

        % compute current Y
        curr_Yopt = 1 / (bbeta * dtempVdinv) + inv_next;
        curr_kopt = (curr_Yopt ./ (((1 - ttau) * (1 - aalpha)) + aalpha) / exp(zz)).^(1/aalpha) * l_ss^((aalpha - 1) / aalpha);

        % compute current value
        curr_g = ttau * (1 - aalpha) * exp(zz) * curr_k^aalpha * l_ss^(1 - aalpha);
        curr_V = (1-bbeta) * (log(curr_Yopt - inv_next) + params.eeta * log(curr_g) - (l_ss^2) / 2) + bbeta * tempV;


        % interpolating the value function given the grid (k, inv) for the shock curr_exst
        Vf_curr = grid.LinearInterpFunction(obj.grid, tempV);



    end


end



function tempV = tempVfunc(inv_next, obj, Vf, curr_st)
    % This is the Vhat variable defined in the main text.

    npts = obj.grid.Npt;

    % unpack some parameters
    ddelta = obj.params.ddelta;
    pphi = obj.params.pphi;

    % current state
    curr_exst = curr_st(1);
    curr_k = curr_st(2);
    curr_inv = curr_st(3);
    pr_next = obj.exogenv.transmat(curr_exst, :);

    k_next = (1 - ddelta) * curr_k + (1 - pphi .* ((inv_next ./ curr_inv) - 1).^2) .* inv_next;

    nextstmat = [repmat((1:obj.exogenv.exnpts), npts, 1), k_next, inv_next];   % [exst_next, k_next, inv_next] all npts by exnpts
    nextstmat = reshape(nextstmat, npts, obj.exogenv.exnpts, 3);
    nextstpts = reshape(nextstmat, npts * obj.exogenv.exnpts, 3);

    V_next = Vf.evaluateAt(nextstpts);
    V_next = reshape(V_next, npts, obj.exogenv.exnpts);

    % compute expected value
    tempV = V_next * pr_next';

end




function V = pv_old(sol, obj, Vf)
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
