% Two examples of how to use the integrate-and-fire simulator. The first
% one contains a single cell and demonstrates how to produce an animation
% of that cell's behavior. The second example defines three cells and
% demonstrates how to define inhibitory and excitatory connections between
% pairs of cells.

%% Example 1: One cell + animation

nsim = nsimulator(1);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);

% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.5, 200, [1]);

%% Example 2: Three cells + connections

nsim = nsimulator(3);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, ie);

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
nsim = nsim.set_mpotential(nsim, 2, -60);

% Specify connections
nsim = nsim.make_excitatory_connection(nsim, 2, 1);
nsim = nsim.make_inhibitory_connection(nsim, 1, 2);

% Simulate and plot
norecordtime = 0.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2 3])
