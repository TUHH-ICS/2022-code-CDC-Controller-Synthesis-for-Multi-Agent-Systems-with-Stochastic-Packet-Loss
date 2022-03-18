clear

addpath('analysis', 'synthesis')
addpath('graphs', 'simulation', 'util')
addpath(genpath('mas-simulation/lib'))

%% Settings
% The script will compare controller synthesized for these probabilites.
% For p = NaN, a nominal controller design according to Massioni will be
% used, for 0 <= p <= 1 a lossy controller.
p_swp = [linspace(0.1, 1, 37), NaN]';

% Simulation parameters
simconf = struct;
simconf.Tf      = 1e5;  % Maximum simulation time [s]
simconf.tol     = 1e-5; % Tolerance before stopping the simulation
simconf.samples = 10;

%% Problem definition
m   = 1;  % Model mass [kg]
b   = 10; % Coefficient of friction [kg/s]
dT  = 1;  % Sampling time [s]
dim = 1;  % Problem dimension

% Discretized state-space model of a mass with friction
A = [ 0   1   ;
      0  -b/m ];
B = [ 0 ; 1/m ];
C = [ 1  0 ];
G = c2d(ss(kron(eye(dim), A), kron(eye(dim), B), kron(eye(dim), C), 0), dT);

% Controller tuning
R = kron(eye(dim), 1);

% Communicaton structure
graph = line_graph(25, 2, false);

% Network config
sinrconf                  = SinrConfiguration();
sinrconf.slotCount        = 10;
sinrconf.wirelessProtocol = WirelessProtocol.underwater_marlin_modem;
sinrconf.pathLoss         = 2.5;
sinrconf.power            = 1e-2;
sinrconf.packetSize       = 6*16;

% Secondary Bernoulli network for reference
bernconf = struct;
bernconf.p     = 0.6;
bernconf.range = Inf;
bernconf.sym   = true;

%% Prepare for synthesis and analysis steps
N = height(graph.Nodes);
simconf.R          = R;   % Penalty on control signal
simconf.dT         = dT;  % Sampling time during simulation
simconf.positions  = [1:N; zeros(dim-1,N)]; % Initial positions of the agents
sinrconf.cycleTime  = dT;  % Time per communication cycle
sinrconf.agentCount = N;

% Assemble the generalized plant
[sysD, sysC, sysP, ny, nu] = prepare_generalized_plant(G, R);

% Calculate graph matrices
L0 = full(laplace_matrix(graph));
A0 = full(adjacency(graph));

%% Try different controllers for the SINR network
H2_sinr   = zeros(size(p_swp));
H2_bern   = zeros(size(p_swp));
H2_nom    = zeros(size(p_swp));
sinr_conv = zeros(size(p_swp), 'logical'); % Convergence flag
bern_conv = zeros(size(p_swp), 'logical'); % Convergence flag

for i = 1:length(p_swp)
    lossy = ~isnan(p_swp(i));
    
    %% Controller synthesis
    tic
    if lossy
        fprintf('Phase %d of %d, p = %g\n', i, length(p_swp), p_swp(i))
        [Kd, Kc] = h2syn_lossy(sysD, sysC, sysP, ny, nu, L0, p_swp(i), 0);
    else
        fprintf('Phase %d of %d, Nominal\n', i, length(p_swp))
        [Kd, Kc] = h2syn_nominal_dual(sysD, addparts(sysC, sysP, 1), ny, nu, L0, 0);
    end
    disp(['Controller synthesis completed in ' format_duration(toc)])

    % Save controller for the Monte-Carlo simulation
    save('controller.mat', 'dT', 'm', 'b', 'Kd', 'Kc')
    
    %% Analyse controller performance without packet loss
    tic
    H2_nom(i) = h2norm_massioni(sysD, addparts(sysC, sysP, 1), Kd, Kc, L0);
    disp(['Nominal analysis completed in ' format_duration(toc)])
    
    %% Monte-Carlo validation        
    tic
    [H2_sinr(i), sinr_conv(i)] = mc_simulate(graph, simconf, sinrconf);
    [H2_bern(i), bern_conv(i)] = mc_simulate(graph, simconf, bernconf);
    disp(['Monte-Carlo simulations completed in ' format_duration(toc)])
end

%% Visualize result
% Filter out results that did not converge
H2_sinr(~sinr_conv) = NaN;
H2_bern(~bern_conv) = NaN;

% Plot results for lossy controllers
figure()
plot(p_swp, H2_sinr, p_swp, H2_bern)
xlabel('Controller Design Transmission Probability')
ylabel('H_2 Performance')
title('H2 Performance under SINR Loss')

legend('Lossy Controller SINR', 'Lossy Controller Bernoulli', 'AutoUpdate', false)

xlim([0, 1])
ylim padded

% Add reference value of nominal controller
hold on
plot(1, H2_nom(isnan(p_swp)), 'x')
H2_sref = H2_sinr(isnan(p_swp));
H2_bref = H2_bern(isnan(p_swp));
if ~isempty(H2_sref)
    plot(1, H2_sref, '+')
end
if ~isempty(H2_bref)
    plot(1, H2_bref, 'd')
end
hold off

%% Export results
name = sprintf('evaluate_sinr_%d.csv', uint32(posixtime(datetime())));
tbl = table(p_swp, H2_sinr, H2_bern, H2_nom);
writetable(tbl, name)