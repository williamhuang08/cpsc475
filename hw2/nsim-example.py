# Two examples of how to use the integrate-and-fire simulator. The first
# one contains a single cell and demonstrates how to produce an animation
# of that cell's behavior. The second example defines three cells and
# demonstrates how to define inhibitory and excitatory connections between
# pairs of cells.

# NOTE: This won't run correctly in a notebook, so you should run the 
# code below (either Example 1 or 2) from a Python prompt.

##### Example 1: One cell + animation
from nsimulator import nsimulator
nsim = nsimulator(1);
ie = 2.5; # nA

# Set initial membrane potential
nsim.set_mpotential(0, -70);
# Define injected current
nsim.set_Ie(0, ie);


# Animate the behavior (display membrane potential vs time)
nsim.animate_integration(0.5, 200, [0])

##### Example 2: Three cells + connections
from nsimulator import nsimulator

nsim = nsimulator(3);
ie = 2.5; # nA

# Define injected current
nsim.set_Ie(0, ie); #setting Ie for neuron 0
nsim.set_Ie(1, ie); #setting Ie for neuron 1

# Set the initial membrane potential
nsim.set_mpotential(0, -60);
nsim.set_mpotential(1, -60);

# Specify connections
nsim.make_excitatory_connection(1, 0);
nsim.make_inhibitory_connection(0, 1);

# Simulate and plot
norecordtime = 0.0; # Integration time (seconds) when not recording 
                    # membrane potentials
recordtime = 0.1; # Integration time for recorded values
nsim.simulate_timeinterval(norecordtime, recordtime);

# Display the recorded values
# 
# Notice that no current is injected into neuron 2 and that neuron is
# completely isolated from the other two (no connections).

nsim.plot_membranepotentials([0, 1, 2])
