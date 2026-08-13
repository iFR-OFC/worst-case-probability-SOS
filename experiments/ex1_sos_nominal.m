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
% time window to explore
T = 1;
% known moments
y1 = 0;
y2 = 0.05;

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

% degree of the hierarchy
list_degree = 1:8;
degree_len = length(list_degree);

% initialize result arrays
results = struct();
results.list_degree = list_degree;
results.built_time  = zeros(degree_len,1);
results.solve_time  = zeros(degree_len,1);
results.upper_bound = zeros(degree_len,1);
results.status      = false(degree_len,1);

for kd =1:degree_len
    %% SETUP OPTIMIZATION PROBLEM FOR EACH ORDER OF THE HIERARCHY
    % get degree
    deg = list_degree(kd);
    
    % decision variables
    w = casos.PS.sym('w', monomials([x;t], 0:deg));
    % sos multipliers
    s  = casos.PS.sym('s', monomials([x;t],0:ceil(deg/2)), 6, 'gram');
    s0 = casos.PS.sym('s0', monomials(x, 0:ceil(deg/2)), 1, 'gram');
    
    % lin and sos decision variables
    x_lin = [w; gamma; v];
    x_sos = [s; s0];
    
    g_lin = [];
    g_sos = [-s(1)*constraint_time-s(2)*g+w-1; 
             -s(3)*constraint_time-s(4)*gX+w;
             -s(5)*constraint_time-s(6)*gX-(nabla(w,t)+nabla(w,x)*fdyn);
             -s0*constraint_X0-subs(w,t,0)-v-gamma(1)*x-gamma(2)*x*x];
    
    % cost function: min -<muk, 1>
    f_cost = gamma(1)*y1+gamma(2)*y2+v;
    
    % Setup the struct
    sos = struct();
    % constraints: linear cone constraints
    sos.g = g_sos;
    % decision variables: lin and sos variables
    sos.x = [x_lin; x_sos];
    % cost function
    sos.f = -f_cost;
    
    % Provide the problem size (size of cones)
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
    % Generate a CaSoS solver instance
    sdp_solver = 'mosek'; 

    tic
    S = casos.sossol('S', sdp_solver, sos, opts); 
    build_time = toc;
    
    %% SOLVE THE OPTIMIZATION AND STORE RESULTS
    tic
    sol = S('lbg', 0, 'ubg', 0);
    solve_time = toc;
    
    % get status of solver
    status = S.stats.UNIFIED_RETURN_STATUS;
    fprintf('d=%d: prob_upper: %f   (%s)\n', deg, full(sol.f), status)
    
    % store results
    results.built_time(kd) = build_time;
    results.solve_time(kd) = solve_time;
    results.upper_bound(kd) = full(sol.f);
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
mean(samples)
var(samples)

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
filename = 'data/ex1_sos_nominal_results.mat';
save(filename, 'results');
