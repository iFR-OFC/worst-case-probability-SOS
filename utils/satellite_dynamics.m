function [X, x0_norm, reni, scale, rk4_step, N_steps, h_norm]=satellite_dynamics(sat_opts)
% Set up the dynamics, RK4, and others

if nargin<1
    sat_opts.plot  = 0;           % 0-no plot / 1-plot
    sat_opts.w     = 3.986e14;    % Earth's gravitational parameter [m^3/s^2]
    sat_opts.scale = [7e6; 3e6; 5e3; 5e3];
    
    % parameters in the real dimensions without normalization
    sat_opts.x0 = [6828e3; 0; 0; 5.4e3];
    sat_opts.t0 = 0;
    sat_opts.tf = 500;
    sat_opts.h  = 1;
    sat_opts.taylor_deg = 2;  
else
    if ~isfield(sat_opts,'plot'),       sat_opts.plot = 0;                       end
    if ~isfield(sat_opts,'w'),          sat_opts.w = 3.986e14;                   end
    if ~isfield(sat_opts,'x0'),         sat_opts.x0 = [6828e3; 0; 0; 5.4e3];     end
    if ~isfield(sat_opts,'scale'),      sat_opts.scale = [7e6; 3e6; 5e3; 5e3];   end
    if ~isfield(sat_opts,'t0'),         sat_opts.t0 = 0;                         end
    if ~isfield(sat_opts,'tf'),         sat_opts.tf = 500;                       end
    if ~isfield(sat_opts,'h'),          sat_opts.h  = 1;                         end
    if ~isfield(sat_opts,'taylor_deg'), sat_opts.taylor_deg = 2;                 end
    
end

%% SYSTEM DEFINITION
% state vector: X = [x; y; vx; vy]
% dynamics: Two-body gravitational problem

% define symbolic variables
x  = casadi.MX.sym('x', 2, 1); % position [x; y]
v  = casadi.MX.sym('v', 2, 1); % velocity [vx; vy]
Xs = [x; v];

% physical parameters
w = sat_opts.w;    % Earth's gravitational parameter [m^3/s^2]

% dynamics: x_dot = v, v_dot = -mu*r/r^3
d = (x'*x)^(3/2);
f = [v(1); v(2); -w*x(1)/d; -w*x(2)/d];

%% SIMULATION PARAMETERS
% integration settings
h  = sat_opts.h;    % time step [s]
t0 = sat_opts.t0;   % initial time [s]
tf = sat_opts.tf;   % final time [s]

% number of integration steps
N_steps = ceil((tf - t0) / h);  

% initial conditions (circular orbit around Earth)
x0 = sat_opts.x0;  % [x; y; vx; vy]

%% POLYNOMIAL APPROXIMATION (CaΣoS)
% approximate dynamics using Taylor polynomial expansion

% create individual component functions
f_components = cell(4, 1);
for i = 1:4
    f_components{i} = casadi.Function(['f', num2str(i)], {Xs}, {f(i)});
end

% define indeterminates for polynomial approximation
X = casos.Indeterminates('X', 4, 1);

% expansion point (initial condition)
x0_exp = x0;

% generate Taylor polynomials for each component
sf_components = cell(4, 1);
for i = 1:4
    sf_components{i} = casos.package.core.Polynomial.from_taylor(...
        f_components{i}, X, x0_exp, sat_opts.taylor_deg);
end

% combine into vector polynomial
sf = [sf_components{1}; sf_components{2}; sf_components{3}; sf_components{4}];

%% NORMALIZATION
% normalize states for better numerical conditioning
% scaling factors: [position_x, position_y, velocity_x, velocity_y]
scale = sat_opts.scale;
inv_scale = 1 ./ scale;

% normalize polynomial dynamics: dx_norm/dt = tf * f(x_norm * scale) * inv_scale
reni = tf*subs(sf, X, X.*scale).*inv_scale;
fcn_reni = reni.to_function;

%% RUNGE-KUTTA INTEGRATION
% implement RK4 integration scheme

% symbolic variables for RK4 step
x_sym = casadi.MX.sym('x_sym',4);   % state
h_sym = casadi.MX.sym('h_sym');     % step size

% RK4 stages
k1 = fcn_reni(x_sym(1), x_sym(2), x_sym(3), x_sym(4));
k2 = fcn_reni(x_sym(1) + h_sym/2 * k1(1), ...
              x_sym(2) + h_sym/2 * k1(2), ...
              x_sym(3) + h_sym/2 * k1(3), ...
              x_sym(4) + h_sym/2 * k1(4));
k3 = fcn_reni(x_sym(1) + h_sym/2 * k2(1), ...
              x_sym(2) + h_sym/2 * k2(2), ...
              x_sym(3) + h_sym/2 * k2(3), ...
              x_sym(4) + h_sym/2 * k2(4));
k4 = fcn_reni(x_sym(1) + h_sym * k3(1), ...
              x_sym(2) + h_sym * k3(2), ...
              x_sym(3) + h_sym * k3(3), ...
              x_sym(4) + h_sym * k3(4));

% RK4 update step
x_next = x_sym + h_sym/6 * (k1 + 2*k2 + 2*k3 + k4);
rk4_step = casadi.Function('rk4_step', {x_sym, h_sym}, {x_next});

%% PERFORM SIMULATION
% normalized time integration
h_norm = h / tf;                         % normalized time step
t_span = linspace(t0, tf, N_steps + 1);  % time vector

% initialize state trajectory
x0_norm = x0.*inv_scale;    % normalized initial condition

% allocate space for trajectory
x_norm = zeros(4, N_steps + 1);
x_norm(:,1) = x0_norm;  

% integrate using RK4
for k = 1:N_steps
    x_norm(:,k+1) = full(rk4_step(x_norm(:,k), h_norm)); 
end

% denormalize results
x_traj = x_norm .* scale;

%% VISUALIZATION
if sat_opts.plot == 1
    % plot results
    fighandle_01 = figure(1);
    
    % orbit plot
    subplot(1, 3, 1);
    plot(x_traj(1, :), x_traj(2, :), 'b-', 'LineWidth', 2);
    hold on;
    plot(x0(1), x0(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    plot(x_traj(1, end), x_traj(2, end), 'gs', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('x [m]');
    ylabel('y [m]');
    title('Orbit Trajectory');
    grid on;
    legend('Trajectory', 'Initial', 'Final', 'Location', 'best');
    
    % position over time
    subplot(1, 3, 2);
    plot(t_span, x_traj(1, :), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(t_span, x_traj(2, :), 'b-', 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel('Position [m]');
    title('Position vs Time');
    grid on;
    legend('x(t)', 'y(t)', 'Location', 'best');
    
    % velocity over time
    subplot(1, 3, 3);
    plot(t_span, x_traj(3, :), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(t_span, x_traj(4, :), 'b-', 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel('Velocity [m/s]');
    title('Velocity vs Time');
    grid on;
    legend('vx(t)', 'vy(t)', 'Location', 'best');
    
    fighandle_01.Position = [1640, 2608, 1038, 283];
end
