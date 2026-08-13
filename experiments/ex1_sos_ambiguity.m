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

%% PROBLEM SETUP
x = casos.Indeterminates('x', 1);
t = casos.Indeterminates('t', 1);

% dynamics
fdyn = 0.5*x*(1-x*x);
fdyn_fcn = fdyn.to_function;

% Parameters
% time window
T = 1;
% moment info
Y = [0; 0.05];
dy = [0.01; 0.001];

% define the set of interest K = {x | g(x) >= 0}
g = 10*(1.*x-0.5)*(1-1.*x);

% define the domain to test
gX = (1.*x+2)*(2-1.*x);

% support of the measures
constraint_X0  = (1.*x+0.5)*(0.5-1.*x);
constraint_time = t*(T-1.*t);

% decision variables
gamma = casos.PS.sym('gamma', 2, 1);
v = casos.PS.sym('v', 1, 1);

% Generate a CaSoS solver instance
sdp_solver = 'mosek'; 

ambiguity = [0; 2; 5; 8; 12];
len_amb = length(ambiguity);

% degree of the hierarchy
list_degree = 1:8;
degree_len = length(list_degree);

% initialize result arrays
results = struct();
results.list_degree     = list_degree;
results.list_ambiguity  = ambiguity;
results.built_time      = zeros(degree_len,len_amb);
results.solve_time      = zeros(degree_len,len_amb);
results.upper_bound     = zeros(degree_len,len_amb);
results.status          = false(degree_len,len_amb);

for ambk = 1:len_amb
    % select ambiguity 'level'
    test_ambiguity = ambiguity(ambk);
    % upper and lower bounds for the optimization
    ubalpha = Y+dy*test_ambiguity; 
    lbalpha = Y-dy*test_ambiguity;
    
    for kd =1:degree_len
        %% SETUP OPTIMIZATION PROBLEM FOR EACH ORDER OF THE HIERARCHY
        deg = list_degree(kd);
        
        % decision variables
        w = casos.PS.sym('w', monomials([x;t], 0:deg));
        
        % sos multipliers
        s  = casos.PS.sym('s', monomials([x;t],0:ceil(deg/2)), 6, 'gram');
        s0 = casos.PS.sym('s0', monomials(x, 0:ceil(deg/2)), 1, 'gram');
        
        % additional terms for the ambiguity components
        ta = casos.PS.sym('ta', 2, 1);
        
        % lin and sos decision variables
        x_lin = [w; gamma; v; ta];
        x_sos = [s; s0];
        
        g_lin = [ta-gamma; ta+gamma];
        
        g_sos = [-s(1)*constraint_time-s(2)*g+w-1; 
                 -s(3)*constraint_time-s(4)*gX+w;
                 -s(5)*constraint_time-s(6)*gX-(nabla(w,t)+nabla(w,x)*fdyn);
                 -s0*constraint_X0-subs(w,t,0)-v-gamma(1)*x-gamma(2)*x*x];
        
        % cost function: min -<muk, 1>
        f_cost = v+0.5*(ubalpha+lbalpha)'*gamma-0.5*(ubalpha-lbalpha)'*ta;
        
        % Setup the struct
        sos = struct();
        % constraints: linear cone constraints
        sos.g = [g_lin; g_sos];
        % decision variables: lin and sos variables
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

        % if true, error returns infeasible
        opts.error_on_fail = false;
        
        %% BUILD SOLVER
        tic
        S = casos.sossol('S', sdp_solver, sos, opts);          
        build_time = toc;
        
        %% SOLVE THE OPTIMIZATION AND STORE RESULTS
        tic
        sol = S('lbx', -inf, 'ubx', +inf, 'lbg', 0, 'ubg', +inf);
        solve_time = toc;
        
        % get status of solver
        status = S.stats.UNIFIED_RETURN_STATUS; 
        fprintf('d=%d: prob_upper: %f   (%s)\n', deg, full(sol.f), status)
        
        % store results
        results.built_time(kd, ambk) = build_time;
        results.solve_time(kd, ambk) = solve_time;
        results.upper_bound(kd, ambk) = full(sol.f);
        results.status(kd, ambk) = strcmp(status, 'SOLVER_RET_SUCCESS');
    
    end

end

%% SAVE RESULTS
% add metadata
results.metadata = struct();
results.metadata.timestamp = datestr(now);

% save to file
filename = 'data/ex1_sos_ambiguity_results.mat';
save(filename, 'results');
