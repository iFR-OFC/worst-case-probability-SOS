%% ------------------------------------------------------------------------
%   
%   Paper: 
%
%   Description: 
%                      
%
%   Needed software: - CasADi 3.6 
%                    - CaΣoS
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3. 
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
r = [0.81; 0.84];
g = +0.005 - (1.*[X(1);X(2)]-r)'*(1.*[X(1);X(2)]-r);

% support of the measures
constraint_X0   = 0.1-(1.*X-x0)'*(1.*X-x0); 
constraint_time = t*(1-1.*t); % restrict the time of interest
constraint_X = 4-X'*X;

% if true, error returns infeasible
opts.error_on_fail = false;

% Generate a CaSoS solver instance
sdp_solver = 'mosek'; 

gamma = casos.PS.sym('gamma', 6, 1);
v = casos.PS.sym('v', 1, 1);

% cost function: min -<muk, 1>
f_cost = v+gamma'*[x0; 0.05+x0(1)*x0(1); 0.01+x0(2)*x0(2)];

% degree of the hierarchy
% set list of degrees 
list_degree = 2:8;
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
    
    % decision variables
    w = casos.PS.sym('w', monomials([X;t], 0:deg));
    % sos multipliers
    s  = casos.PS.sym('s', monomials([X;t],0:ceil(deg/2)), 6, 'gram');
    s0 = casos.PS.sym('s0', monomials(X, 0:ceil(deg/2)), 1, 'gram');
    
    % lin and sos decision variables
    x_lin = [w; gamma; v];
    x_sos = [s; s0];
    
    % linear constraints
    g_lin = [];

    % sos constraints
    g_sos = [-s(1)*constraint_time-s(2)*g+w-1;
             -s(3)*constraint_time-s(4)*constraint_X+w;
             -s(5)*constraint_time-s(6)*constraint_X-(nabla(w,t)+nabla(w,X)*reni);
             -s0*constraint_X0-subs(w,t,0)-v-gamma'*[X; X(1)*X(1); X(2)*X(2)]
        ];
    
    % Setup the struct
    sos = struct();
    % constraints: linear cone constraints
    sos.g = [g_lin; g_sos];
    % decision variables: measure variables
    sos.x = [x_lin; x_sos];
    % cost function
    sos.f = -f_cost;
    
    % Provide the problem size i.e. size of cones
    nx_lin = length(x_lin);
    nx_sos = length(x_sos);
    ng_lin = length(g_lin);
    ng_sos = length(g_sos);
    
    % decision variable cones
    opts.Kx.lin = nx_lin;   % sos decision variables
    opts.Kx.sos = nx_sos;   % sos decision variables
    % constraint cones
    opts.Kc.lin = ng_lin;   % linear constraints
    opts.Kc.sos = ng_sos;   % nonnegative measure cone constraint
    
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
    results.upper_bound(kd) = full(sol.f);
    results.status(kd) = strcmp(status, 'SOLVER_RET_SUCCESS');

end

%% SAVE RESULTS
% add metadata
results.metadata = struct();
results.metadata.timestamp = datestr(now);

% save to file
filename = 'data/ex2_sos_nominal_results.mat';
save(filename, 'results');
