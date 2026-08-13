%% ------------------------------------------------------------------------
%   
%   Paper: "Worst-Case Probability Bounds for Finite-Horizon Safety under 
%          Moment Uncertainty"
%
%   Description: Find the worst case probability given certain moments 
%                associated with initial measure using a moment 
%                program. This example considers a unsafe set that is time 
%                varying.
%                      
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% ------------------------------------------------------------------------

% add path to utility folder with auxiliary functions
addpath("utils/")

% set up options for satellite dynamics generation (OPTIONAL)
sat_opts = [];
sat_opts.plot = 0;

% generate scaled satellite dynamics (polynomial approximation)
[X, x0, reni, scale, rk4_step, N, h_norm] = satellite_dynamics(sat_opts); 

%% PROBLEM SETUP
t = casos.Indeterminates('t', 1);

% define the set of interest K = {x | g(x) >= 0}
r0 = [0.81; 0.84]; 
v0 = 10*[-0.01; -0.01];
g = +0.005 - (1.*[X(1);X(2)]-r0-v0*t)'*(1.*[X(1);X(2)]-r0-v0*t);

% time window
T = 1;
% support of the measures
supp_mu_init  = 0.1-(1.*X-x0)'*(1.*X-x0);
supp_mu_term  = t*(T-1.*t);
supp_mu_occu  = t*(T-1.*t);

% SDP solver
sdp_solver = 'mosek';

% degree of the hierarchy
list_degree = 2:5;
degree_len = length(list_degree);

% initialize result arrays
results = struct();
results.list_degree = list_degree;
results.built_time  = zeros(degree_len,1);
results.solve_time  = zeros(degree_len,1);
results.upper_bound = zeros(degree_len,1);
results.status      = false(degree_len,1);

for kd = 1:degree_len
    %% SETUP OPTIMIZATION PROBLEM FOR EACH ORDER OF THE HIERARCHY
    deg = list_degree(kd);
    
    % get degree for mu_init and mu_term
    degv = 2*deg-reni.maxdeg;
    degv = degv + rem(degv,2);
    
    % measure to be on the set K
    muk = casos.PS.sym('muk', monomials([t; X(1:2)],0:degv/2), 'dual');    % measure to get the probability
    mu_init = casos.PS.sym('mu_init', monomials(X,0:degv/2), 'dual');      % initial measure
    mu_occu = casos.PS.sym('mu_occu', monomials([t; X],0:degv/2), 'dual'); % occupation measure
    mu_term = casos.PS.sym('mu_term', monomials([t; X],0:degv/2), 'dual'); % terminal measure
    
    % measure decision variables (moments)
    x_meas = [mu_init.primalize; mu_term.primalize; mu_occu.primalize; muk.primalize];
    
    % set the liouville equation
    mu_liouv = mu_occu.liouville(reni,X,t);
    v = monomials([t; X], 0:mu_liouv.maxdeg);
    liouville_eq = project(mu_term.primalize, v)-project(mu_init.primalize, v)-project(primalize(mu_liouv),v);
    
    g_lin = [liouville_eq;              % dynamics constraint
             mu_init.evaluate(1)-1;         % ensure it is a probability measure
             mu_term.evaluate(1)-1;
             mu_init.evaluate(X(1))-x0(1);  % set the mean of the initial measure (x-axis)
             mu_init.evaluate(X(2))-x0(2);  % set the mean of the initial measure (y-axis)
             mu_init.evaluate(X(3))-x0(3);
             mu_init.evaluate(X(4))-x0(4);
             mu_init.evaluate((1.*X(1)-x0(1))'*(1.*X(1)-x0(1))) - 0.05;
             mu_init.evaluate((1.*X(2)-x0(2))'*(1.*X(2)-x0(2))) - 0.01;
             ];
    
    pmuk = muk.primalize;

    g_meas = [primalize(mu_init.support(supp_mu_init));        % set support of the initial measure    
              primalize(mu_term.support(supp_mu_term));        % set support of the terminal measure
              primalize(mu_occu.support(supp_mu_occu));        % set support of the occupation measure
              primalize(mu_occu.support(4-X'*X)); % ensure some bounded support        
              primalize(muk.support(g));                       % set support of muk
              project(mu_term.primalize,pmuk.sparsity)-pmuk;
              ];
    
    % cost function: min -<muk, 1>
    f_cost = -muk.evaluate(1);
    
    % Setup the struct
    sos = struct();
    % constraints: linear cone constraints
    sos.g = [g_lin; g_meas];
    % decision variables: measure variables
    sos.x = x_meas;
    % cost function
    sos.f = f_cost;
    
    % Provide the problem size i.e. size of cones
    nx_meas = length(x_meas);
    ng_lin = length(g_lin);
    ng_meas = length(g_meas);
    
    % decision variable cones
    opts.Kx.meas = nx_meas;   % measure decision variables
    % constraint cones
    opts.Kc.lin  = ng_lin;    % linear constraints
    opts.Kc.meas = ng_meas;   % nonnegative measure cone constraint

    % if true, error returns infeasible
    opts.error_on_fail = false;
    
    %% BUILD SOLVER
    tic
    S = casos.sossol('S', sdp_solver, sos, opts); 
    build_time = toc;
    
    %% SOLVE THE OPTIMIZATION AND STORE RESULTS
    tic
    sol = S('lbg', 0, 'ubg', 0);
    solve_time = toc;

    % get status of solver
    status = S.stats.UNIFIED_RETURN_STATUS;
    fprintf('d=%d: prob_upper: %f   (%s)\n', deg, -full(sol.f), status)

    % store results
    results.built_time(kd) = build_time;
    results.solve_time(kd) = solve_time;
    results.upper_bound(kd) = -full(sol.f);
    results.status(kd) = strcmp(status, 'SOLVER_RET_SUCCESS');

end

%% SAVE RESULTS
% add metadata
results.metadata = struct();
results.metadata.timestamp = datestr(now);

% save to file
filename = 'data/ex2_meas_moving_results.mat';
save(filename, 'results');