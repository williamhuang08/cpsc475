import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class nsimulator:
    def __init__(self, nneurons):
        """NSIMULATOR An integrate-and-fire neuron simulator

        SIM = nsimulator(NNEURONS) returns an instance of the neuron
        simulator object. The simulation is as described in Chapter 5 of "Theoretical
        Neuroscience" by Peter Dayan and L.F. Abbott, 1st ed. The model follows
        Eq. 5.35 and Eq. 5.43. Where needed, variables defined as in Fig. 5.20.

        NNEURONS is the number of neurons in the simulation.

        SIM is the simulator. It is an object with methods (callable like 
        regular functions). These functions are described below. 
        Note that many of these functions when called will change the state 
        of the simulator object.

        sim.set_Ie(sim, ni, current)  

        Sets the injected current into cell 'ni' to the value of 'current'
        (nA).


        sim.set_mpotential(ni, mpotential)

        Sets the membrane potential of cell 'ni' to 'mpotential' (mV).
        Initially the membrane potentials are initialized to a random value.


        sim.make_excitatory_connection(i, j)

        Defines a synaptic connection with presynaptic cell 'i' and
        postsynaptic cell 'j', i.e. 'i' affects 'j'.


        sim.make_inhibitory_connection(i, j)

        Defines a synaptic connection with presynaptic cell 'i' and
        postsynaptic cell 'j', i.e. 'i' affects 'j'.


        sim.simulate_timeinterval(norecordtime, recordtime)

        Simulates behavior first for a time given by 'norecordtime' and then
        for a time given by 'recordtime' during which a record is kept of the
        membrane potentials for each neuron. 'recordtime' must be positive
        and nonzero (time in seconds).


        sim.plot_membranepotentials(neurons)

        Plots the membrane potentials recorded during a simulation. Function
        must be called after a simulate_timeinterval() call. 'neurons' is a
        vector of neuron indices that should be plotted. 
        E.g., sim.plot_membranepotentials([0 2 3]) will plot the records
        for neurons 0, 2, and 3 in plots appearing in a new figure.

        sim.animate_integration(recordtime, frames, neurons)

        Same as the plotting function above but this one performs the
        simulation and displays the results 'frames' times.

        January, 2008.
        Written by Pavel Dimitrov for CPSC475 at Yale University.
        """
        self.nneurons = nneurons
        self.Vs = -20e-3*np.random.random(nneurons) - 50e-3 # initial membrane potentials
        self.Ies = np.zeros((nneurons,1)) # injected current
        
        self.PsM = np.zeros((nneurons, nneurons))    # Matrix of Ps values PsM(i,j) is the Ps in cell j incoming from cell i
        self.PstimeM = 0.1*np.ones((nneurons, nneurons)) # Matrix of time values for the calculation of Ps(t) (eq. 5.35)
        self.EsM = np.zeros((nneurons, nneurons))     # Matrix of synapse current (0mV excitatory and -80mV is inhibitory)
        self.dt = 1e-4

    def set_Ie(self, ni, current):
        # ni  -- neuron number
        # current -- new injected current value (nA)
    
        self.Ies[ni] = current*1e-9 

    def set_mpotential(self, ni, mpotential):
        # ni  -- neuron number
        # mpotential -- new membrane potential in mV
    
        self.Vs[ni] = mpotential*1e-3 
    
    def make_excitatory_connection(self, i,j):
        if i < 0 or j < 0 or i >= self.nneurons or j >= self.nneurons:
            raise ValueError('i and j must be valid neuron indices, i.e. between 0 and %d-1', self.nneurons)
    
        self.PsM[i,j] = 1
        self.EsM[i,j] = 0  # Volts
    

    def make_inhibitory_connection(self, i, j):
        if i < 0 or j < 0 or i >= self.nneurons or j >= self.nneurons:
            raise ValueError('i and j must be valid neuron indices, i.e. between 0 and %d-1', self.nneurons)
    
        self.PsM[i,j] = 1
        self.EsM[i,j] = -80*1e-3 # Volts

    def simulate_timestep(self):
        dt = self.dt
        
        Pmax = 5
        tau_s = 10e-3 # time scale for Ps
        
        Vreset = -80*1e-3    # V
        Vth =    -54*1e-3    # V
        Vmax =    10*1e-3    # V

        E_L = -70*1e-3 # (V) leak voltage
        Rm  =  10*1e6  # ohms, total membrane resistance
        tau =  20*1e-3 # s, tau=cm*rm % time constant
        
        self.Vmax = Vmax
        
        for ni in range(self.nneurons):
            V =  self.Vs[ni]
            Ie = self.Ies[ni]

            if V == Vmax: # is this an action potential
                V = Vreset
                
                # action potential -> time-since-AP=0
                self.PstimeM[ni, self.PsM[ni,:] != 0] = 0
            else:
                synapsecontribution = 0
                if np.flatnonzero(self.PsM).size > 0:
                    for nbn in np.flatnonzero(self.PsM[:,ni]):
                        Ps = self.PsM[nbn,ni]
                        Es = self.EsM[nbn,ni]
                        synapsecontribution -=  0.05*Ps*(V-Es) 

                
                dVdt =  1./tau*( E_L - V + Rm*Ie + synapsecontribution) 
                V = V + dt*dVdt

                if V > Vth:
                    V = Vmax

            
            self.Vs[ni] = V

        
        # now update the Ps and the time-since-action-potential PstimeM
        mask = (self.PsM != 0).astype('double')
        self.PstimeM += dt*mask
        tM = self.PstimeM
        self.PsM = mask * (Pmax*tM/tau_s*np.exp(1. - tM/tau_s))
        

    def simulate_timeinterval(self, etint, tint):
        dt = self.dt
        
        if tint < 0:
           raise ValueError('recordtime must be nonzero and positive')
        
        if etint > 0:
            timesteps = int(round(etint/dt))
            for ti in range(timesteps):
                self.simulate_timestep()

        
        timesteps = int(round(tint/dt))
        t = np.zeros(timesteps+1)
        VrecordM = np.zeros((self.nneurons, timesteps))
        
        for ti in range(timesteps):
            # first record the current membrane potentials for each neuron
            for ni in range(self.nneurons):
                VrecordM[ni,ti] = self.Vs[ni]


            self.simulate_timestep()
            t[ti+1] = t[ti] + dt


        self.VrecordM = VrecordM.copy()
        self.timev = t[:-1].copy()


    def animate_integration(self, tint, frames, neuronlist):
        dt = self.dt

        nlsz = len(neuronlist)

        timesteps = int(round(tint/dt))
        t = np.zeros(timesteps+1)
        VrecordM = np.zeros((self.nneurons, timesteps))

        tstp100 = int(round(timesteps/frames))
        # simulate
        for ti in range(timesteps):
            # first record the current membrane potentials for each neuron
            for ni in range(self.nneurons):
                VrecordM[ni,ti] = self.Vs[ni]
            
            self.simulate_timestep()
            t[ti+1] = t[ti] + dt

        self.VrecordM = VrecordM.copy()
        self.timev = t[:-1].copy()

        fig, axes = plt.subplots(nlsz,1,figsize=(6,2.5*nlsz))
        fig.set_tight_layout(True)
        fig.subplots_adjust(hspace=.5)
        if nlsz == 1:
            axes = [axes]
        trains = []

        for dispi,ax in enumerate(axes):
            # initialize plot
            ax.set_xlabel('time (ms)')
            ax.set_xlim(0,tint*1000)

            ax.set_ylabel('membrane potential (mV)')
            v = VrecordM[neuronlist[dispi], :]
            ax.set_ylim(v.min()*1000,v.max()*1000)
            ax.tick_params(labelsize=8)
            
            
            ax.set_title('Neuron {0} ({1} spikes)'.format(neuronlist[dispi],
                np.sum(v==self.Vmax)))


            trains.append(ax.plot([], [])[0])



        def plot_mps(i):

            t = self.timev[:i+1]
            for dispi in range(nlsz):
                v = self.VrecordM[neuronlist[dispi], :i+1]
                assert v.size >= 1
                trains[dispi].set_data(t*1000, v*1000)

            return trains

        ani = FuncAnimation(fig, plot_mps, frames=range(0,timesteps+1,tstp100), interval=100, repeat=False)
        plt.show()

    def plot_membranepotentials(self, neuronlist):
        nlsz = len(neuronlist)
        tint = (len(self.timev)-1)*self.dt
        fig, axes = plt.subplots(nlsz,1,figsize=(5,2.5*nlsz))
        if nlsz == 1:
            axes = [axes]
        else:
            fig.set_tight_layout(True)
            fig.subplots_adjust(hspace=.5)
        for dispi in range(nlsz):
            ax = axes[dispi]
            ax.set_xlabel('time (ms)')
            ax.set_ylabel('membrane potential (mV)')

            t = self.timev
            v = self.VrecordM[neuronlist[dispi], :]
            ax.plot(t*1000,v*1000)
            ax.set_title('Neuron {0} ({1} spikes)'.format(neuronlist[dispi],
                np.sum(v==self.Vmax)))
        plt.show()



