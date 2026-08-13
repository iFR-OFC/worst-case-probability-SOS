%% ------------------------------------------------------------------------
%   
%   Paper: "Worst-Case Probability Bounds for Finite-Horizon Safety under 
%          Moment Uncertainty"
%
%   Description: Find the worst case probability given certain moments on 
%                the initial measure using a moment program.
%                     
%   License: GNU GENERAL PUBLIC LICENSE Version 3
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
% known moments     
y1 = 0;
y2 = 0.05;

% define the set of interest K = {x | g(x) >= 0}
g = (1.*x-0.5)*(1-1.*x);

% define the domain to test
gX = (1.*x+2)*(2-1.*x);

% support of the measures
supp_mu_init  = (1.*x+0.5)*(0.5-1.*x);
supp_mu_term  = t*(T-1.*t); 
supp_mu_occu  = t*(T-1.*t); 

% degree of the hierarchy
list_degree = 3:8;
degree_len = length(list_degree);

% SDP solver
sdp_solver = 'mosek'; 

% initialize result arrays
results = struct();
results.list_degree = list_degree;
results.built_time  = zeros(degree_len,1);
results.solve_time  = zeros(degree_len,1);
results.upper_bound = zeros(degree_len,1);
results.status      = false(degree_len,1);

for kd = 1:degree_len
    %% SETUP OPTIMIZATION PROBLEM FOR EACH ORDER OF THE HIERARCHY
    % get degree
    deg = list_degree(kd);
    
    % get degree for mu_init and mu_term
    degv = 2*deg-fdyn.maxdeg;
    degv = degv+rem(degv,2);
    
    % measure to be on the set K
    muk     = casos.PS.sym('muk', monomials(x,0:degv/2), 'dual');          % auxiliary measure 
    mu_init = casos.PS.sym('mu_init', monomials(x,0:degv/2), 'dual');      % initial measure
    mu_occu = casos.PS.sym('mu_occu', monomials([t; x],0:degv/2), 'dual'); % occupation measure
    mu_term = casos.PS.sym('mu_term', monomials([t; x],0:degv/2), 'dual'); % terminal measure
    
    % measure decision variables
    x_meas = [mu_init.primalize; 
              mu_term.primalize; 
              mu_occu.primalize; 
              muk.primalize];
    
    % Liouville equation
    mu_liouv = mu_occu.liouville(fdyn,x,t);
    v = monomials([t;x], 0:mu_liouv.maxdeg);
    liouville_eq = project(mu_term.primalize, v)-...
                   project(mu_init.primalize, v)-...
                   project(primalize(mu_liouv),v);
   
    % linear constraints
    g_lin = [liouville_eq;                  % dynamics constraint
             mu_init.evaluate(1)-1;         % ensure mu_init is a probability measure
             mu_term.evaluate(1)-1;         % ensure mu_term is a probability measure
             mu_init.evaluate(x)-y1;        % set the mean of the initial measure (x-axis)
             mu_init.evaluate(x*x)-y2;      % set noncentral 2nd order moment
             mu_term.dot(t)-1;              % terminal meaure should on t=1
             ];
    
    % primalize muk
    pmuk = muk.primalize;
    
    % measure constraints
    g_meas = [primalize(mu_init.support(supp_mu_init));        % set support of the initial measure    
              primalize(mu_term.support(supp_mu_term));        % set support of the terminal measure
              primalize(mu_term.support(gX));                  % set support of the terminal measure
              primalize(mu_occu.support(supp_mu_occu));        % set support of the occupation measure
              primalize(mu_occu.support(gX));                  % ensure some bounded support        
              primalize(muk.support(g));                       % set support of muk
              project(mu_term.primalize,pmuk.sparsity)-pmuk;   % domination constraint
              ];
    
    % cost function: min -<muk, 1>
    f_cost = -muk.evaluate(1);
    
    % setup the struct
    sos = struct();
    % constraints: linear cone constraints
    sos.g = [g_lin; g_meas];
    % decision variables: measure variables
    sos.x = x_meas;
    % cost function
    sos.f = f_cost;
    
    % Provide the problem size (size of cones)
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

%% GET A LOWER BOUND WITH MONTE CARLO SIMULATIONS   
% polynomial density function: p(x) = -6*x^2+0*x+1.5
p = [-6, 0, 1.5]; 

% define the bounds
x_min = -0.5;
x_max = 0.5;

% generate a fine grid of points to integrate the PDF numerically
x_grid = linspace(x_min, x_max, 10000);
pdf_values = polyval(p, x_grid);
% normalize the PDF so it integrates to 1 over [x_min, x_max]
area = trapz(x_grid, pdf_values);
pdf_normalized = pdf_values / area;
% compute the Cumulative Distribution Function (CDF)
cdf_values = cumtrapz(x_grid, pdf_normalized);

% generate Uniformly distributed random numbers for sampling (0 to 1)
num_samples = 1e4;
U = rand(num_samples, 1);

% use interpolation to find the x-values corresponding to U (Inverse Transform)
samples = interp1(cdf_values, x_grid, U, 'pchip');

% verify the moments, namely the first and second moment as a sanity check
fprintf('Empirical mean: %f\nEmpirical variance: %f\n', mean(samples), var(samples))

% run simulation on each sample
dt = 0.001;
N = T/dt;

Xtest = zeros(length(samples), N+1);
Xtest(:,1) = samples;
prob = zeros(1, N);
for i=1:N
    Xtest(:,i+1) = Xtest(:,i) + dt*full(fdyn_fcn(Xtest(:,i)));
    % get the probability from Xtest at the current instance
    idx1 = Xtest(:,i+1)>=0.5;
    idx2 = Xtest(:,i+1)<=1;
    prob(i) = sum(idx1&idx2)/length(samples);
end
results.lower_bound = prob;

%% SAVE RESULTS
% add metadata
results.metadata = struct();
results.metadata.timestamp = datestr(now);

% save to file
filename = 'data/ex1_meas_nominal_results.mat';
save(filename, 'results');
