% Author: Jiacheng Li
% Date: Oct. 20, 2024
%
% This is the main script for HW 1 of the course ECON 8210.



% ---------------------------------------------------
% 0. Housekeeping
% ---------------------------------------------------
close all; clear all;
delete(gcp('nocreate'));
parpool(4)

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


% Initial guess for the endogenous variables [k, l, c]
x0 = [1, 0.5, 0.5]';

% Solve the steady state equations using fsolve
options = optimoptions('fsolve', 'Display', 'iter');
[x_ss, fval, exitflag] = fsolve(@(x) SSeq(x, params), x0, options);
[~, inv_ss, g_ss, y_ss, V_ss] = SSeq(x_ss, params);

% Display the results
k_ss = x_ss(1);
l_ss = x_ss(2);
c_ss = x_ss(3);

disp('Steady state values:');
disp(['Capital (k): ', num2str(k_ss)]);
disp(['Labor (l): ', num2str(l_ss)]);
disp(['Consumption (c): ', num2str(c_ss)]);
disp(['Investment (inv): ', num2str(inv_ss)]);
disp(['Output (y): ', num2str(y_ss)]);
disp(['Gov consumption (g): ', num2str(g_ss)]);
disp(['Value: ', num2str(V_ss)]);

ssvals = struct();
ssvals.k_ss = k_ss;
ssvals.l_ss = l_ss;
ssvals.c_ss = c_ss;
ssvals.inv_ss = inv_ss;
ssvals.y_ss = y_ss;
ssvals.g_ss = g_ss;
ssvals.V_ss = V_ss;








% ---------------------------------------------------
% 2. VFI with fixed grid
% ---------------------------------------------------


% --------------- prepare ----------------
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

% transmat_ttau = 1;      % a deterministic version to test grid bound
% transmat_zz = 1;
% ttau_grid = 0.25;
% zz_grid = 0;


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
exogenv.exgrid = 1:exogenv.exnpt;


% everything into obj
obj = struct('params', params, 'exogenv', exogenv, 'ssvals', ssvals);





% --------------- VFI fixed grid ----------------

disp('------------------------------------');
disp('VFI with fixed grid');
disp('------------------------------------');

% set up grid points
ggrid_start = struct('kNpt', 250, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_start, ssvals, exogenv);
obj.nskip = 1;  % first, do not apply acceleration step

% initial guess
Vmat0 = V_ss * ones(obj.ggrid.maingrid.Npt, 1);
Pmat0 =  [c_ss, l_ss, inv_ss, k_ss, g_ss] .* ones(obj.ggrid.maingrid.Npt, 1);


% run VFI
tic;
[Vf1, Pf1] = runVFI(obj, Vmat0, Pmat0);
toc;

% plot result
plotres(obj, Vf1, Pf1, 'VFI with fixed grid');


% print EE error
EEerr(obj, Vf1, Pf1);





% ---------------------------------------------------
% 3. VFI with switching between policy and VF iteration
% ---------------------------------------------------

obj.nskip = 10;  % now apply acceleration step

% run VFI
tic;
[Vf1, Pf1] = runVFI(obj, Vmat0, Pmat0);
toc;

% plot result
plotres(obj, Vf1, Pf1, 'VFI with acceleration');









% ---------------------------------------------------
% 4. VFI with muiltigrid grid
% ---------------------------------------------------

% --------------- VFI start round ----------------
% set up grid points
ggrid_start = struct('kNpt', 100, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_start, ssvals, exogenv);


% initial guess
Vmat0 = V_ss * ones(obj.ggrid.maingrid.Npt, 1);
Pmat0 =  [c_ss, l_ss, inv_ss, k_ss, g_ss] .* ones(obj.ggrid.maingrid.Npt, 1);


% run VFI
tic;
[Vf1, Pf1] = runVFI(obj, Vmat0, Pmat0);
toc;


% --------------- VFI fine round ----------------
% set up grid points
ggrid_fine = struct('kNpt', 500, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_fine, ssvals, exogenv);

% initial guess
Vmat0 = Vf1.evaluateAt(obj.ggrid.maingrid.Pointmat)';
Pmat0 = Pf1.evaluateAt(obj.ggrid.maingrid.Pointmat)';


% run VFI
tic;
[Vf2, Pf2] = runVFI(obj, Vmat0, Pmat0);
toc;


% --------------- VFI final round ----------------
% set up grid points
ggrid_finest = struct('kNpt', 5000, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_finest, ssvals, exogenv);

% initial guess
Vmat0 = Vf2.evaluateAt(obj.ggrid.maingrid.Pointmat)';
Pmat0 = Pf2.evaluateAt(obj.ggrid.maingrid.Pointmat)';


% run VFI
tic;
[Vf3, Pf3] = runVFI(obj, Vmat0, Pmat0);
toc;


% plot result
plotres(obj, Vf3, Pf3, 'VFI with multigrid (final round)');


% print EE error
EEerr(obj, Vf3, Pf3);


% --------------- plot IRF -----------------
% positive tax shock
shock_id = [3, 3];
transpath1 = simul_IRF(shock_id, 30, obj, Vf3, Pf3);

% positive productivity shock
shock_id = [2, 5];
transpath2 = simul_IRF(shock_id, 30, obj, Vf3, Pf3);













% ---------------------------------------------------
% 5. VFI with stochastic grid
% ---------------------------------------------------
% set up grid points
ggrid_stoch = struct('kNpt', 500, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_stoch, ssvals, exogenv);

% initial guess
Vmat0 = V_ss * ones(obj.ggrid.maingrid.Npt, 1);
Pmat0 =  [c_ss, l_ss, inv_ss, k_ss, g_ss] .* ones(obj.ggrid.maingrid.Npt, 1);
obj.nskip = 10;

% run VFI
tic;
[Vf4, Pf4] = runstogVFI(obj, Vmat0, Pmat0);
toc;


% plot result
plotres(obj, Vf4, Pf4, 'VFI with stochastic grid');


% print EE error
EEerr(obj, Vf4, Pf4);






% ---------------------------------------------------
% 6. Endogenous grid method
% ---------------------------------------------------
% -------- small VFI to get some initial guess with slope ----------
% set up grid points
ggrid_start = struct('kNpt', 250, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_start, ssvals, exogenv);
obj.nskip = 10;

% initial guess
Vmat0 = V_ss * ones(obj.ggrid.maingrid.Npt, 1);
Pmat0 =  [c_ss, l_ss, inv_ss, k_ss, g_ss] .* ones(obj.ggrid.maingrid.Npt, 1);

% run VFI
tic;
[Vf1, Pf1] = runVFI(obj, Vmat0, Pmat0);
toc;


% -------- perform EGM with fixed labor --------
% set up grid points
ggrid_stoch = struct('kNpt', 500, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_stoch, ssvals, exogenv);

% initial guess
Vmat1 = Vf1.evaluateAt(obj.ggrid.maingrid.Pointmat)';
Pmat1 = Pf1.evaluateAt(obj.ggrid.maingrid.Pointmat)';

% run VFI
tic;
[Vf2, Pf2] = runVFI_EGM(obj, Vmat1, Pmat1);
toc;


% -------------- feed to standard VFI ------------
% set up grid points
ggrid_finest = struct('kNpt', 5000, 'invNpt', 50);
obj.ggrid = constructgrid(ggrid_finest, ssvals, exogenv);
obj.nskip = 10;

% initial guess
Vmat3 = Vf2.evaluateAt(obj.ggrid.maingrid.Pointmat)';
Pmat3 = Pf2.evaluateAt(obj.ggrid.maingrid.Pointmat)';


% run VFI
tic;
[Vf3, Pf3] = runVFI(obj, Vmat3, Pmat3);
toc;


plotres(obj, Vf3, Pf3);








% ---------------------------------------------------
% Useful functions
% ---------------------------------------------------

function [res, inv, g, y, V] = SSeq(x, params)

    % Unpack the parameters
    aalpha  = params.aalpha;    % capital share
    bbeta   = params.bbeta;     % discount factor
    ddelta  = params.ddelta;    % depreciation rate
    eeta    = params.eeta;    % gov consumption weight

    ttauSS  = params.ttauSS;    % steady state tax on labor
    zzSS    = params.zzSS;      % steady state productivity
    
    ppsi = (1 - ttauSS) * (1 - aalpha) + aalpha;
    

    % Unpack the endogenous variables
    k = x(1);   % capital
    l = x(2);   % labor
    c = x(3);   % consumption

    
    % intermediate variables
    y = exp(zzSS) * k^aalpha * l^(1-aalpha);    % output
    
    % Compute the steady state residuals
    res = ones(3, 1);

    res(1) = l - (1 / c) * ppsi * (1 - aalpha) * y/l;
    res(2) = ppsi * aalpha * y / k - (1/bbeta - 1 + ddelta);
    res(3) = c + ddelta * k - ppsi * y;

    inv = ddelta * k;
    g = ttauSS * (1 - aalpha) * y;
    V = log(c) + eeta * log(g) - (l^2) / 2;
end



function ggrid = constructgrid(ggrid, ssvals, exogenv)

    % get steady state values
    k_ss = ssvals.k_ss;
    inv_ss = ssvals.inv_ss;

    ggrid.kgrid     = linspace(0.1 * k_ss, 2 * k_ss, ggrid.kNpt);
    ggrid.invgrid   = linspace(0.3 * inv_ss, 1.7 * inv_ss, ggrid.invNpt);

    ggrid.maingrid = grid.TensorGrid({exogenv.exgrid, ggrid.kgrid, ggrid.invgrid});
    
end




function [Vf, Pf] = runVFI(obj, Vmat0, Pmat0)

    % hyperparameters
    maxIter = 1000;
    tol = 1e-6;
    Verr = 100;

    Vmat = Vmat0;
    Pmat = Pmat0;

    iter = 0;


    % value function iteration
    while iter < maxIter && Verr > tol

        if mod(iter, obj.nskip) == 0
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

    % interpolate the value function
    Vf = grid.LinearInterpFunction(obj.ggrid.maingrid, Vmat);
    Pf = grid.LinearInterpFunction(obj.ggrid.maingrid, Pmat);

end



function [Vf1, Pf1] = runstogVFI(obj, Vmat0, Pmat0)

    % benchmark grid
    benchgrid = obj.ggrid.maingrid;

    % interpolate the value function
    Vf0 = grid.LinearInterpFunction(benchgrid, Vmat0);
    Pf0 = grid.LinearInterpFunction(benchgrid, Pmat0);

    % hyperparameters
    maxIter = 1000;
    tol = 1e-6;
    Verr = 100;

    iter = 0;


    % value function iteration
    while iter < maxIter && Verr > tol

        if mod(iter, obj.nskip) == 0
            accel = 0;
        else
            accel = 1;
        end

        % uniformly draw sample k grid
        kgridsample = sort(randsample(obj.ggrid.kgrid, 200));
        obj.ggrid.maingrid = grid.TensorGrid({obj.exogenv.exgrid, kgridsample, obj.ggrid.invgrid});

        % reevaluate Vmat, Pmat
        Vmat = Vf0.evaluateAt(obj.ggrid.maingrid.Pointmat)';
        Pmat = Pf0.evaluateAt(obj.ggrid.maingrid.Pointmat)';
        
        % main update
        [Vmat, Pmat] = updateV(Vmat, Pmat, obj, accel);

        % reinterpolate
        Vf1 = grid.LinearInterpFunction(obj.ggrid.maingrid, Vmat);
        Pf1 = grid.LinearInterpFunction(obj.ggrid.maingrid, Pmat);

        % evaluate at the benchmark to compute the distance
        Vmat = Vf1.evaluateAt(benchgrid.Pointmat)';

        Verr = max(abs(Vmat - Vmat0));
        disp(['Iteration ', num2str(iter), ', Vf error: ', num2str(Verr), ', Accelerate: ', num2str(accel)]);
        
        Vf0 = Vf1;
        Pf0 = Pf1;
        Vmat0 = Vmat;
        iter = iter + 1;
    end

    % reset the grid
    obj.ggrid.maingrid = benchgrid;
end


function [] = EEerr(obj, Vf, Pf)

    bbeta = obj.params.bbeta;
    pphi = obj.params.pphi;

    % unpack grid points
    gridmat = obj.ggrid.maingrid.Pointmat;

    % get the policies
    polmat = Pf.evaluateAt(gridmat);
    cgrid = polmat(1, :)';
    invpgrid = polmat(3, :)';
    kpgrid = polmat(4, :)';

    Err = zeros(obj.ggrid.maingrid.Npt, 1);

    % get the expected value
    for ii = 1:obj.ggrid.maingrid.Npt

        state = gridmat(ii, :);

        % next state
        state_next = [(1:obj.exogenv.exnpt)', [kpgrid(ii), invpgrid(ii)] .* ones(obj.exogenv.exnpt, 1)];

        % approximate the derivative
        h = 1e-6;
        state_next_invp = [(1:obj.exogenv.exnpt)', [kpgrid(ii), invpgrid(ii) + h] .* ones(obj.exogenv.exnpt, 1)];
        state_next_kp = [(1:obj.exogenv.exnpt)', [kpgrid(ii) + h, invpgrid(ii)] .* ones(obj.exogenv.exnpt, 1)];

        Vf_next = Vf.evaluateAt(state_next);
        Vf_next_invp = Vf.evaluateAt(state_next_invp);
        Vf_next_kp = Vf.evaluateAt(state_next_kp);
        Vf_next_invp = (Vf_next_invp - Vf_next) / h;
        Vf_next_kp = (Vf_next_kp - Vf_next) / h;

        EVf_next_invp = obj.exogenv.transmat(state(1), :) * Vf_next_invp';
        EVf_next_kp = obj.exogenv.transmat(state(1), :) * Vf_next_kp';

        dkpdinvp = 1 - pphi * (invpgrid(ii) / state(3) - 1).^2 + 2 * pphi * (invpgrid(ii) / state(3)) .* (invpgrid(ii) / state(3) - 1);

        c = cgrid(ii);
        Err(ii) = 1 - bbeta * (EVf_next_invp + EVf_next_kp * dkpdinvp) / ((1-bbeta) * 1/c);
    end

    disp('====================================');
    % Display the percentiles of the Euler equation errors
    percentiles = prctile(Err, [0, 25, 50, 75, 100]);
    disp('Percentiles of the Euler equation errors:');
    disp(['Min: ', num2str(percentiles(1))]);
    disp(['25th percentile: ', num2str(percentiles(2))]);
    disp(['50th percentile: ', num2str(percentiles(3))]);
    disp(['75th percentile: ', num2str(percentiles(4))]);
    disp(['Max: ', num2str(percentiles(5))]);
    meanErr = mean(Err);
    disp(['Mean Euler equation error: ', num2str(meanErr)]);
    disp('====================================');
end



function [Vf, Pf] = runVFI_EGM(obj, Vmat0, Pmat0)

    % hyperparameters
    maxIter = 1000;
    tol = 1e-6;
    Verr = 100;

    Vmat = Vmat0;
    Pmat = Pmat0;

    iter = 0;


    % value function iteration
    while iter < maxIter && Verr > tol
        
        [Vmat, Pmat] = updateV_EGM(Vmat, Pmat, obj);

        Verr = max(abs(Vmat - Vmat0));
        disp(['Iteration ', num2str(iter), ', Vf error: ', num2str(Verr)]);
        
        Vmat0 = Vmat;
        iter = iter + 1;
    end

    % interpolate the value function
    Vf = grid.LinearInterpFunction(obj.ggrid.maingrid, Vmat);
    Pf = grid.LinearInterpFunction(obj.ggrid.maingrid, Pmat);

end



function [Vmat, Pmat] = updateV(Vmat, Pmat, obj, accel)
    % Each iteration in fixed grid VFI

    % Unpack the parameters
    aalpha = obj.params.aalpha;

    % unpack grid points
    gridmat = obj.ggrid.maingrid.Pointmat;
    Npt = obj.ggrid.maingrid.Npt;
    exstmat = obj.exogenv.exstmat;

    % interpolate the value function
    Vf = grid.LinearInterpFunction(obj.ggrid.maingrid, Vmat);

    parfor ii = 1:Npt

        % local variables
        state = gridmat(ii, :);


        xopt = Pmat(ii, :);

        if ~accel
            % unpack the state
            ttau = exstmat(state(1), 1);
            zz = exstmat(state(1), 2);
            k = state(2);
    
            % upper bound on consumption
            cmax = ((1 - ttau) * (1 - aalpha) + aalpha) * exp(zz) * k^aalpha / (1 - aalpha)^((aalpha - 1) / 2);
            obj2min = @(c) - pVc(c, gridmat(ii, :), Vf, obj);
    
            % minimize - value
            opts = optimset('Display', 'off');
            xopt(1) = fminbnd(obj2min, 0.3, cmax, opts);
        end

        [pV, xopt(2), xopt(3), xopt(4), xopt(5)] = pVc(xopt(1), gridmat(ii, :), Vf, obj);

        Vmat(ii) = pV
        Pmat(ii, :) = xopt;
    end

end



function [val, lopt, inv_next, k_next, g] = pVc(c, state, Vf, obj)

    % unpack params and states
    aalpha = obj.params.aalpha;
    ddelta = obj.params.ddelta;
    pphi = obj.params.pphi;
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






function [Vmat, Pmat] = updateV_EGM(Vmat, Pmat, obj)

    % Unpack the parameters
    aalpha = obj.params.aalpha;
    bbeta = obj.params.bbeta;
    l_ss = obj.ssvals.l_ss;

    % unpack grid points
    Npt = obj.ggrid.maingrid.Npt;
    exstmat = obj.exogenv.exstmat;
    invpgrid = obj.ggrid.invgrid;     % use for policy
    invpNpt = length(invpgrid);
    
    % store the total extended endogenous grid points for interpolation later
    endogrid = zeros(invpNpt, 3, Npt);
    Vendogrid = zeros(invpNpt, Npt);
    Pendogrid = zeros(invpNpt, 5, Npt);

    % interpolate the value function
    Vf = grid.LinearInterpFunction(obj.ggrid.maingrid, Vmat);

    for ii = 1:Npt

        % rowids = (ii-1)*invpNpt+1:ii*invpNpt;

        % local variables
        state = obj.ggrid.maingrid.Pointmat(ii, :);
        inv = state(3);
        tau = exstmat(state(1), 1);
        zz = exstmat(state(1), 2);
        ppsi = (1 - tau) * (1 - aalpha) + aalpha;

        % evaluate the derivative using finite difference
        [Vhatgrid, kpgrid] = Vhat(invpgrid, state, Vf, obj);

        h = 0.00001;
        dVhatdinvpgrid = (Vhat(invpgrid + h, state, Vf, obj) - Vhatgrid) / h;

        Ygrid = 1 ./ (bbeta./(1-bbeta) .* dVhatdinvpgrid) + invpgrid;
        kgrid = (Ygrid ./ (ppsi * exp(zz) * l_ss^(1 - aalpha))).^(1 / aalpha);
        cgrid = Ygrid - invpgrid;
        g_grid = tau * (1 - aalpha) * exp(zz) * kgrid.^aalpha * l_ss^(1 - aalpha);

        pVgrid = (1-bbeta) * (log(cgrid) - l_ss^2 / 2) + bbeta * Vhatgrid;


        % endogenous grid points
        endogrid(:, :, ii) = [state(1) * ones(invpNpt, 1), kgrid', inv * ones(invpNpt, 1)];
        Vendogrid(:, ii) = pVgrid;

        % policy
        Pendogrid(:, :, ii) = [cgrid', l_ss * ones(invpNpt, 1), invpgrid', kpgrid', g_grid'];
    end

    endogrid = reshape(permute(endogrid, [1,3,2]), invpNpt*Npt, 3);
    Vendogrid = reshape(permute(Vendogrid, [1,3,2]), invpNpt*Npt, 1);
    Pendogrid = reshape(permute(Pendogrid, [1,3,2]), invpNpt*Npt, 5);

    % interpolate the value/policy function
    Vf = scatteredInterpolant(endogrid, Vendogrid, 'linear', 'boundary');
    Vmat = Vf(obj.ggrid.maingrid.Pointmat);

    for jj = 1:5
        Pf = scatteredInterpolant(endogrid, Pendogrid(:, jj), 'linear', 'boundary');
        Pmat(:, jj) = Pf(obj.ggrid.maingrid.Pointmat);
    end
end




function [Vhatgrid, kpgrid] = Vhat(invpgrid, state, Vf, obj)

    % Unpack the parameters
    aalpha = obj.params.aalpha;
    pphi = obj.params.pphi;
    ddelta = obj.params.ddelta;

    % unpack states
    k = state(2);
    inv = state(3);
    tau = obj.exogenv.exstmat(state(1), 1);
    zz = obj.exogenv.exstmat(state(1), 2);
    ppsi = (1 - tau) * (1 - aalpha) + aalpha;

    invpNpt = length(invpgrid);
    exnpt = obj.exogenv.exnpt;
    exgrid = obj.exogenv.exgrid;

    % kp grid
    kpgrid = (1 - ddelta) * k + (1 - pphi * (invpgrid ./ inv - 1).^2) .* invpgrid;

    % next period value and expectation
    nextstate3d = repmat([kpgrid', invpgrid'], 1, 1, exnpt);
    nextstate3d = [repmat(reshape((1:exnpt), 1, 1, exnpt), invpNpt, 1), nextstate3d];
    nextstate = reshape(permute(nextstate3d, [1, 3, 2]), exnpt * invpNpt, 3);

    Vpgrid = reshape(Vf.evaluateAt(nextstate), exnpt, invpNpt);
    Vhatgrid = obj.exogenv.transmat(state(1), :) * Vpgrid;
end






function plotres(obj, Vf, Pf, name)

    ggrid = obj.ggrid;
    exogenv = obj.exogenv;

    % read some steady state values
    inv_ss = obj.ssvals.inv_ss;
    k_ss = obj.ssvals.k_ss;
    c_ss = obj.ssvals.c_ss;
    g_ss = obj.ssvals.g_ss;
    l_ss = obj.ssvals.l_ss;


    grid_plot = [exogenv.meanstid * ones(length(ggrid.kgrid), 1), ...
        ggrid.kgrid', inv_ss * ones(ggrid.kNpt, 1)];
    

    % Plot the policy functions
    pol_plot = Pf.evaluateAt(grid_plot);
    
    f = figure('visible','off');
    subplot(3, 2, 1);
    plot(ggrid.kgrid, Vf.evaluateAt(grid_plot), 'LineWidth', 2);
    title('Value Function');
    xlabel('Capital (k)');
    ylabel('Value Function (V)');
    grid on;
    
    subplot(3, 2, 2);
    plot(ggrid.kgrid, pol_plot(1, :), 'LineWidth', 2);
    hold on;
    plot(k_ss, c_ss, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Policy Function for Consumption');
    xlabel('Capital (k)');
    ylabel('Consumption (c)');
    grid on;
    
    subplot(3, 2, 3);
    plot(ggrid.kgrid, pol_plot(2, :), 'LineWidth', 2);
    hold on;
    plot(k_ss, l_ss, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Policy Function for Labor');
    xlabel('Capital (k)');
    ylabel('Labor (l)');
    grid on;
    
    subplot(3, 2, 4);
    plot(ggrid.kgrid, pol_plot(3, :), 'LineWidth', 2);
    hold on;
    plot(k_ss, inv_ss, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Policy Function for Investment');
    xlabel('Capital (k)');
    ylabel('Investment (inv)');
    grid on;

    subplot(3, 2, 5);
    plot(ggrid.kgrid, pol_plot(4, :), 'LineWidth', 2);
    hold on;
    plot(k_ss, k_ss, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    plot(ggrid.kgrid, ggrid.kgrid, 'k--', 'LineWidth', 2);
    hold off;
    title('Policy Function for Next period Capital');
    xlabel('Capital (k)');
    ylabel('Capital (kp)');
    grid on;
    
    subplot(3, 2, 6);
    plot(ggrid.kgrid, pol_plot(5, :), 'LineWidth', 2);
    hold on;
    plot(k_ss, g_ss, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Policy Function for Government spending');
    xlabel('Capital (k)');
    ylabel('Gov spending (g)');
    grid on;

    sgtitle(name);

    % Save the plot to a PNG file
    saveas(f, [name, '.png']);
    
end 






function transpath = simul_IRF(shock_id, NT, obj, Vf, Pf)
    % Simulate the impulse response function
    % INPUTs
    %     sol_ss: a 1D array of steady state values [k, l, c, inv]
    %     shock_id: a 1D array of shock id [ttau_id, zz_id]
    %     NT: a scalar of number of periods to simulate
    %     obj: a structure containing the model parameters, grids and exogenous environment
    % OUTPUTs
    %    transpath: a 2D array of simulated path [exstid, k, inv, c, l, V, ttau, zz]

    % Unpack the steady state values
    k_ss = obj.ssvals.k_ss;
    l_ss = obj.ssvals.l_ss;
    c_ss = obj.ssvals.c_ss;
    inv_ss = obj.ssvals.inv_ss;
    g_ss = obj.ssvals.g_ss;
    V_ss = obj.ssvals.V_ss;


    % get shock
    shock_st = find(obj.exogenv.indlist(:, 1) == shock_id(1) & obj.exogenv.indlist(:, 2) == shock_id(2));

    % simulate the path
    transpath = zeros(NT+1, 9);   % [exstid, k, inv, c, l, V, ttau, zz]
    transpath(:, 1) = obj.exogenv.meanstid;
    transpath(:, 2) = k_ss;
    transpath(:, 3) = inv_ss;
    transpath(:, 4) = c_ss;
    transpath(:, 5) = l_ss;
    transpath(:, 6) = g_ss;
    transpath(:, 7) = V_ss;    % value function

    transpath(:, 8) = obj.params.ttauSS;
    transpath(:, 9) = obj.params.zzSS;


    for tt = 1:NT
        if tt == 2  % initial state
            transpath(tt, 1) = shock_st;
            transpath(tt, 8) = obj.exogenv.exstmat(shock_st, 1);
            transpath(tt, 9) = obj.exogenv.exstmat(shock_st, 2);
        end

        curr_st = transpath(tt, 1:3);

        % policy path
        pol_sol = Pf.evaluateAt(curr_st);
        transpath(tt, 4) = pol_sol(1);
        transpath(tt, 5) = pol_sol(2);
        transpath(tt, 6) = pol_sol(5);
        transpath(tt, 7) = Vf.evaluateAt(curr_st);

        % state transition
        transpath(tt + 1, 2:3) = [pol_sol(4), pol_sol(3)];
    end

    transpath = transpath(1:NT, :);

    % plot transpath in subplots
    f = figure('visible','off');
    subplot(3, 3, 1);
    plot(transpath(:, 1), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 1), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Exogenous State Index');
    xlabel('Time');
    ylabel('Index');
    grid on;

    subplot(3, 3, 2);
    plot(transpath(:, 2), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Capital (k)');
    xlabel('Time');
    ylabel('k');
    grid on;

    subplot(3, 3, 3);
    plot(transpath(:, 3), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 3), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Investment (inv)');
    xlabel('Time');
    ylabel('inv');
    grid on;

    subplot(3, 3, 4);
    plot(transpath(:, 4), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 4), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Consumption (c)');
    xlabel('Time');
    ylabel('c');
    grid on;

    subplot(3, 3, 5);
    plot(transpath(:, 5), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 5), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Labor (l)');
    xlabel('Time');
    ylabel('l');
    grid on;

    subplot(3, 3, 6);
    plot(transpath(:, 6), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 6), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Government Spending (g)');
    xlabel('Time');
    ylabel('g');
    grid on;

    subplot(3, 3, 7);
    plot(transpath(:, 7), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 7), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Value Function (V)');
    xlabel('Time');
    ylabel('V');
    grid on;

    subplot(3, 3, 8);
    plot(transpath(:, 8), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 8), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Tax Rate (ttau)');
    xlabel('Time');
    ylabel('ttau');
    grid on;

    subplot(3, 3, 9);
    plot(transpath(:, 9), 'LineWidth', 2);
    hold on;
    plot(1, transpath(1, 9), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;
    title('Productivity Shock (zz)');
    xlabel('Time');
    ylabel('zz');
    grid on;

    name = ['IRF_shock_', mat2str(shock_id)];

    saveas(f, [name, '.png']);

    
end