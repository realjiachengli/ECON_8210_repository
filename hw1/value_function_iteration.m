% Author: Jiacheng Li
% Date: Oct. 20, 2024
%
% This is the main script for HW 1 of the course ECON 8210. The Miranda-Fackler compecon toolbox is used
% to solve some of the problems.



% ---------------------------------------------------
% 0. Housekeeping
% ---------------------------------------------------
close all; clear all;

addpath('../compecon/CEtools');



% ---------------------------------------------------
% 1. Compute the steady state
% ---------------------------------------------------

% Define the parameters
params.aalpha = 0.33;
params.bbeta = 0.97;
params.ddelta = 0.9;
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
params.transmat_ttau = [0.9, 0.1, 0; 0.05, 0.9, 0.05; 0, 0.1, 0.9];
params.transmat_zz = [
        0.9727 0.0273 0      0      0;
        0.0041 0.9806 0.0153 0      0;
        0      0.0082 0.9836 0.0082 0;
        0      0      0.0153 0.9806 0.0041;
        0      0      0      0.0273 0.9727
    ];
ttau_grid = [0.2, 0.25, 0.3];
zz_grid = [-0.0673, -0.0336, 0, 0.0336, 0.0673];
params.transmat = kron(params.transmat_ttau, params.transmat_zz);

% define grid
k_grid = linspace(0.7 * k_ss, 1.3 * k_ss, 250);
inv_grid = linspace(inv_ss - 0.5 * k_ss, inv_ss + 0.5 * k_ss, 50);

[KK, II, ZZ, TT] = ndgrid(k_grid, inv_grid, zz_grid, ttau_grid);
params.pointmat = [KK(:), II(:), ZZ(:), TT(:)];





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

    res(1) = (1 - ddelta) * k - inv;
    res(2) = eeta * (1 - aalpha) * (1 / l) - (l - (1 / c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * (1 - aalpha) * y);
    res(3) = c + inv - ((1 - ttauSS) * (1 - aalpha) + aalpha) * y;
    res(4) = (1 / c) - bbeta * (eeta * aalpha * (1 / k) + (1/c) * ((1 - ttauSS) * (1 - aalpha) + aalpha) * aalpha * y);


end


function transpath = simul_IRF(Vf, Pf, params)
end


function [Vf, Pf] = solve_VFI_fixed_grid(sol_ss, params)
    % Solve the model using value function iteration with fixed grid
    % INPUTs
    %     params: a structure containing the model parameters
    % OUTPUTs
    %     Vf: a 1D array of value functions [V]
    %     Pf: a 2D array of policy functions [c, l]

    % unpack the grids
    pointmat = params.pointmat;
    transmat = params.transmat;

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
    Vf0 = ones(size(pointmat, 1), 1) * ...
        (log(c_ss) + eeta * log(ttauSS * (1 - aalpha) * exp(zzSS) * k_ss^aalpha * l_ss^(1 - aalpha)) - (l_ss^2) / 2);

    

end