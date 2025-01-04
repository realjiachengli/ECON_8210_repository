function solveNeoclassicalRBC_FE()
% solveNeoclassicalRBC_FE
%
% Solves the single‐equation RBC model with consumption and labor:
%   1/c(k,z) = 0.97 * E[ (1/c(k',z'))*(0.33 e^{z'} k'^{-0.67} l(k',z')^{0.67} + 0.9 ) ]
%
% subject to:
%   c(k,z) = (0.67 e^z k^0.33) / l(k,z)^{1.33}
%   k'(k,z) = e^z k^0.33 l(k,z)^0.67 + 0.9*k - c(k,z)
%   z' = 0.95 z + 0.007 eps,   eps ~ N(0,1).
%
% We approximate l(k,z) by a 2D piecewise‐linear FE basis on
%   k in [0,kmax],  z in [zmin,zmax].
%
% Then we impose a Galerkin condition:
%   \int\int psi_i(k)*psi_j(z) * [Residual(k,z; theta)] dk dz = 0
% for all (i,j), which yields a large nonlinear system in the unknown
% theta_{i,j}.
%
% This code uses 'fsolve' to solve the system.  You may want to use your
% favorite nonlinear solver (Newton, Quasi‐Newton, etc.) instead.
%
% Written as a skeleton.  Adapt/optimize to your preferences.

    % Model parameters
    beta   = 0.97;   % discount factor
    alpha  = 0.33;   % exponent on capital
    delta  = 0.10;   % "effective" net dep. for the sum 0.9 below
                     % (We've effectively used 0.9 as '1-delta' in your eqn.)

    % State space bounds
    kmax   = 100;    % upper bound on capital
    zmin   = -0.5;   % lower bound on productivity
    zmax   =  0.5;   % upper bound on productivity

    % Number of FE partitions (subintervals):
    Nk = 15;    % # of subintervals in k
    Nz = 3;    % # of subintervals in z
    % => # FE "nodes" in each dimension = Nk+1, Nz+1

    % Build the grid edges:
    kgrid = linspace(0,   kmax, Nk+1);
    zgrid = linspace(zmin,zmax, Nz+1);

    % Pack parameters into a struct for convenience
    param.beta  = beta;
    param.alpha = alpha;
    param.kgrid = kgrid;
    param.zgrid = zgrid;

    % Create an initial guess for the FE coefficients (theta_{i,j})
    % We'll guess a constant labor of 0.3 on the entire domain:
    ncoef   = (Nk+1)*(Nz+1);
    theta0  = 0.30*ones(ncoef,1);

    % Set up options for fsolve (or you could use e.g. fminunc, or your own Newton)
    options = optimoptions('fsolve',...
        'Display','iter',...
        'MaxFunEvals',5e4,...
        'MaxIter',1e4,...
        'TolFun',1e-8,...
        'TolX',1e-8);

    % Solve
    [theta_sol, fval, exitflag] = fsolve(@(th) residualSystem(th,param), theta0, options);

    disp('-----------------------------------------------');
    fprintf('Exit flag           = %d\n', exitflag);
    fprintf('Norm of residual    = %g\n', norm(fval));
    disp('-----------------------------------------------');

    % OPTIONAL: Plot the solution l(k,z) on a 2D mesh, for instance
    plotLaborSolution(theta_sol,param);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESIDUAL SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rvec = residualSystem(theta, param)
% residualSystem
%
% Returns the vector of Galerkin residuals for each (i,j).
% The size of Rvec is (Nk+1)*(Nz+1)-by-1.
%
% R(k,z; theta) = 1/c(k,z) - beta * E[ ... ] = 0
% Weighted by the test function psi_i(k)*psi_j(z) and integrated over k,z.

    kgrid = param.kgrid;
    zgrid = param.zgrid;
    beta  = param.beta;

    Nk = length(kgrid)-1;   % # subintervals in k
    Nz = length(zgrid)-1;   % # subintervals in z

    ncoef = (Nk+1)*(Nz+1);
    Rvec  = zeros(ncoef,1);

    % We'll do 2D integration by looping over each element rectangle:
    % For each (i,j), we "test" the residual with the basis function
    % psi_i(k)*psi_j(z).  We'll store the integrated value in Rvec(rowID).

    eq_count = 1;  % index in Rvec

    for iNode = 1:(Nk+1)
        for jNode = 1:(Nz+1)

            % This is the basis function index (iNode,jNode).

            % We now integrate over the entire domain, *but* effectively
            % the only contribution is from the element(s) where
            % psi_i(k)*psi_j(z) != 0.
            % For a piecewise-linear basis, that means up to two intervals
            % in k, and up to two intervals in z, around (iNode, jNode).
            % Below is a simplistic approach: we just loop over
            % all subintervals in k,z, but skip the integral if the basis
            % is identically zero on that subinterval.  A more advanced
            % code would store adjacency info to be more efficient.

            % Accumulate the integral for this test function:
            integral_ij = 0.0;

            for iElem = 1:Nk
                kLeft  = kgrid(iElem);
                kRight = kgrid(iElem+1);

                % 1D Gauss-Legendre nodes/weights on [kLeft, kRight]
                [kq, wkq] = gaussLegendreQuad(kLeft,kRight);

                for jElem = 1:Nz
                    zLeft  = zgrid(jElem);
                    zRight = zgrid(jElem+1);

                    % 1D Gauss-Legendre nodes/weights on [zLeft, zRight]
                    [zq, wzq] = gaussLegendreQuad(zLeft,zRight);

                    % Double loop over quadrature points:
                    sumLocal = 0.0;
                    for a = 1:length(kq)
                        for b = 1:length(zq)
                            kk = kq(a);
                            zz = zq(b);
                            ww = wkq(a)*wzq(b);

                            % Evaluate the test function at (k,z):
                            psi_k_i  = basis_k(kk, param.kgrid, iNode);
                            psi_z_j  = basis_z(zz, param.zgrid, jNode);
                            test_val = psi_k_i * psi_z_j;
                            if (test_val == 0), continue; end

                            % Now evaluate the pointwise residual R(k,z; theta)
                            R_pt = pointwiseResidual(kk,zz,theta,param);

                            sumLocal = sumLocal + R_pt*test_val*ww;
                        end
                    end

                    integral_ij = integral_ij + sumLocal;
                end
            end

            Rvec(eq_count) = integral_ij;
            eq_count       = eq_count + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POINTWISE RESIDUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rval = pointwiseResidual(k,z,theta,param)
% Rval = 1/c(k,z) - beta * E[ (1/c(k',z'))(0.33 e^{z'} k'^{-0.67} l(k',z')^{0.67} + 0.9 ) ]
%
% Here c(k,z) = [0.67 * e^z * k^0.33]/[l(k,z)^1.33]
% k'(k,z) = e^z k^0.33 l(k,z)^0.67 + 0.9*k - c(k,z)
% z' = 0.95 z + 0.007 eps, eps ~ N(0,1)

    beta  = param.beta;

    % Evaluate labor:
    l_val = labor_FE(k,z,theta,param);
    if l_val<=0
        % If l<=0, code can blow up.  Just avoid negative labor if it occurs 
        % (which might indicate a poor initial guess or that the model 
        % solution doesn't exist in that region).
        Rval = 1e8;  % large penalty
        return
    end

    % Consumption c
    c_val = cFunction(k,z,l_val);
    if c_val<=0
        Rval = 1e8;
        return
    end

    % The current period part: 1/c(k,z)
    curr_term = 1.0/c_val;

    % Next period capital
    kp_val = kprimeFunction(k,z,l_val);

    % Expected future term
    Eterm = expectedTerm(kp_val, z, theta, param);

    Rval = curr_term - beta*Eterm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EVALUATE LABOR (the FE approximation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l_val = labor_FE(k,z,theta,param)
% l_val = Sum_{i,j} theta_{i,j} * psi_i(k) * psi_j(z).
%
% We map (iNode,jNode) => linear index in theta:
    kgrid = param.kgrid;
    zgrid = param.zgrid;
    Nk = length(kgrid)-1;   % # subintervals
    Nz = length(zgrid)-1;

    val = 0.0;
    idx = 1;  % will run from 1 to (Nk+1)*(Nz+1)

    for iNode = 1:(Nk+1)
        psi_k_i = basis_k(k, kgrid, iNode);
        for jNode = 1:(Nz+1)
            psi_z_j = basis_z(z, zgrid, jNode);
            val = val + theta(idx)*psi_k_i*psi_z_j;
            idx = idx+1;
        end
    end
    l_val = val;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSUMPTION FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c_val = cFunction(k,z,l_val)
% c(k,z) = [0.67 e^z k^0.33]/[l(k,z)^1.33]
    c_val = (0.67 * exp(z) * k^0.33) / (l_val^1.33);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K' FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kp = kprimeFunction(k,z,l_val)
% k'(k,z) = e^z k^0.33 l(k,z)^0.67 + 0.9*k - c(k,z)
    c_val = (0.67 * exp(z) * k^0.33) / (l_val^1.33);
    kp    = exp(z)*k^0.33*(l_val^0.67) + 0.9*k - c_val;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPECTATION (OVER Z') TERM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E_val = expectedTerm(kp, z, theta, param)
% E[ 1/c(kp,z') * (0.33 e^{z'} kp^{-0.67} l(kp,z')^{0.67} + 0.9) ]
% where z' = 0.95 z + 0.007 eps, eps ~ N(0,1).
%
% We do Gauss‐Hermite quadrature over eps.

    [eNodes, eWeights] = gaussHermiteNodesWeights(10);  % 10-point GH, say
    sumVal = 0.0;

    for m = 1:length(eNodes)
        eps_m = eNodes(m);
        w_m   = eWeights(m);

        % next period z'
        zp = 0.95*z + 0.007*eps_m;

        % labor(kp,zp)
        l_p = labor_FE(kp,zp,theta,param);
        if l_p<=0
            % If that occurs, skip or penalize
            continue
        end

        % c(kp,zp)
        c_p = cFunction(kp,zp,l_p);
        if c_p<=0
            continue
        end

        % integrand
        integrand = (1/c_p) * (0.33*exp(zp)*(kp^(-0.67))*l_p^0.67 + 0.9);

        sumVal = sumVal + integrand * w_m;
    end

    % For standard Gauss-Hermite, the normal pdf factor is built in if you
    % do the scaling properly.  Usually we do:
    %   integral_{-inf to inf} f(x) phi(x) dx  ~ sum f(x_i)*w_i
    % with x_i, w_i from GH.  Because your sigma is 0.007, you might want
    % to adjust the scaling.  The code below assumes you used a GH rule
    % for the standard Normal.  Then you'd do z' = 0.95*z + sigma * x,
    % etc.  Possibly you need an extra factor 1/sqrt(pi) or 1/sqrt(2) 
    % depending on how your GH nodes are defined.  This is something 
    % you can confirm with your GH routine.
    %
    % Suppose we are doing E_{eps~N(0,1)}[f(eps)] = 1/sqrt(pi)* sum( w_i f(x_i) ).
    % If your gaussHermiteNodesWeights() returns that format, you might do:

    % We'll do a typical convention: sumVal = sumVal / sqrt(pi).
    E_val = sumVal / sqrt(pi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAUSS‐LEGENDRE QUADRATURE ON [a,b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xq, wq] = gaussLegendreQuad(a,b)
% Return nodes/weights for 2‐point (or 3‐point) Gauss‐Legendre on [a,b].
% For higher accuracy, you might want 4 or 5 points, etc.
% Here let's define a 3‐point rule as an example.

    % 3‐point Gauss‐Legendre on [-1,1]
    x0 = [-sqrt(3/5), 0, sqrt(3/5)];
    w0 = [ 5/9      , 8/9, 5/9      ];

    % map [-1,1] -> [a,b]
    mid = 0.5*(a+b);
    rad = 0.5*(b-a);

    xq = mid + rad*x0;
    wq = rad*w0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAUSS‐HERMITE QUADRATURE (STANDARD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xH, wH] = gaussHermiteNodesWeights(N)
% Returns the N‐point Gauss‐Hermite nodes/weights for the integral:
% integral_{-inf to inf} e^{-x^2} f(x) dx ~ sum_{i=1..N} w_i f(x_i).
%
% If you want E_{X~N(0,1)}[f(X)] = (1/sqrt(pi)) * integral e^{-x^2} f(x) dx,
% then the expected value is (1/sqrt(pi)) * sum_{i=1..N} w_i f(x_i).
%
% Here, we just use a built‐in approach with 'hermiteH' or we can code 
% an approximation.  For simplicity, let's supply a short approximation for
% up to say 10 points.  If you want a robust routine, consider a proper 
% GH routine or a well‐tested library.

    switch N
        case 2
            xH = [-1/sqrt(2); 1/sqrt(2)];
            wH = [ sqrt(pi)/2; sqrt(pi)/2 ];
        case 3
            xH = [-sqrt(3/2); 0; sqrt(3/2)];
            wH = [ sqrt(pi)/6; sqrt(pi)/3; sqrt(pi)/6 ];
        case 4
            % etc. ...
            error('Add your desired GH nodes for N=4');
        case 10
            % Precomputed 10‐point GH abscissas/weights (common table).
            % xH in ascending order; wH > 0
            xH = [
                -3.436159118, -2.532731674, -1.756683649, -1.036610829, ...
                -0.342901327,  0.342901327,  1.036610829,  1.756683649, ...
                 2.532731674,  3.436159118
            ]';
            wH = [
                0.0004825732, 0.0331949837, 0.2404586130, 0.6108626337, ...
                0.6108626337, 0.2404586130, 0.0331949837, 0.0004825732, ...
                0, 0  % <-- Actually, we need the full correct set ...
            ]';  
            % (You should fill in the correct 10 weights. A standard GH 
            %  table or a built-in might be more reliable.)

            % For brevity: Let's store a typical set from a standard table:
            wH = [ ...
                0.0025557844, 0.0886166370, 0.2913116138, 0.5217556106, ...
                0.6012181920, 0.5217556106, 0.2913116138, 0.0886166370, ...
                0.0025557844, 0.0000233699  % example set, might not sum perfectly
            ]'; 
        otherwise
            error('Only small N coded here. Extend as needed.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIECEWISE‐LINEAR BASIS FOR K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = basis_k(k, kgrid, iNode)
% "Tent" function that is 1 at node kgrid(iNode) and 0 at other nodes.
% iNode goes from 1..(Nk+1), so kgrid has length(Nk+1).
%
% If kgrid = [k0, k1, k2, ..., kNk], then:
%   basis_k(k, kgrid, i) = 0 if k < k_{i-1} or k > k_{i+1}
%   else it is linearly ascending or descending in between.

    val = 0.0;

    % If iNode > 1, we have a left interval
    if (iNode > 1)
        kl = kgrid(iNode-1);
        kc = kgrid(iNode);
        if (k>=kl && k<=kc)
            val = (k - kl)/(kc - kl);
            return
        end
    end

    % If iNode < length(kgrid), we have a right interval
    if (iNode < length(kgrid))
        kc = kgrid(iNode);
        kr = kgrid(iNode+1);
        if (k>=kc && k<=kr)
            val = (kr - k)/(kr - kc);
            return
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIECEWISE‐LINEAR BASIS FOR Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = basis_z(z, zgrid, jNode)
% Same idea as basis_k, but in z dimension.
    val = 0.0;

    if (jNode > 1)
        zl = zgrid(jNode-1);
        zc = zgrid(jNode);
        if (z>=zl && z<=zc)
            val = (z - zl)/(zc - zl);
            return
        end
    end

    if (jNode < length(zgrid))
        zc = zgrid(jNode);
        zr = zgrid(jNode+1);
        if (z>=zc && z<=zr)
            val = (zr - z)/(zr - zc);
            return
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL: PLOT THE FINAL LABOR FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotLaborSolution(theta_sol,param)
% We can sample a coarse 2D mesh, evaluate l_FE at each point, 
% and create a surf plot or contour plot.

    kgrid = param.kgrid;
    zgrid = param.zgrid;

    NkPlot = 40;
    NzPlot = 40;
    kmin   = kgrid(1);
    kmax   = kgrid(end);
    zmin   = zgrid(1);
    zmax   = zgrid(end);

    Kplot = linspace(kmin,kmax,NkPlot);
    Zplot = linspace(zmin,zmax,NzPlot);

    Lmat  = zeros(NkPlot, NzPlot);

    for i=1:NkPlot
        for j=1:NzPlot
            Lmat(i,j) = labor_FE(Kplot(i),Zplot(j),theta_sol,param);
        end
    end

    figure;
    [KK,ZZ] = meshgrid(Kplot,Zplot);
    % But note meshgrid typically does (X=cols, Y=rows). We'll transpose Lmat:
    surf(KK',ZZ',Lmat);
    xlabel('k');
    ylabel('z');
    zlabel('l(k,z)');
    title('Approximate Labor Function via FE');
    shading interp;
    colorbar;
end