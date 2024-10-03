function sim = nsimulator(nneurons)
%NSIMULATOR An integrate-and-fire neuron simulator
%   SIM = nsimulator(NNEURONS) returns an instance of the neuron
%   simulator. The simulation is as described in Chapter 5 of "Theoretical
%   Neuroscience" by Peter Dayan and L.F. Abbott, 1st ed. The model follows
%   Eq. 5.35 and Eq. 5.43. Where needed, variables defined as in Fig. 5.20.
%
%   NNEURONS is the number of neurons in the simulation.
%
%   SIM is the simulator. It is a structure with fields containing
%   function handles (callable like regular functions). These functions
%   are described below. Note that many of these functions return a
%   modified simulator structure. This must be assigned to your simulator
%   variable for the changes to take effect.
%
%   sim = sim.set_Ie(sim, ni, current)  
%   
%     Sets the injected current into cell 'ni' to the value of 'current'
%     (nA).
%   
%   
%   sim = sim.set_mpotential(sim, ni, mpotential)
%   
%     Sets the membrane potential of cell 'ni' to 'mpotential' (mV).
%     Initially the membrane potentials are initialized to a random value.
%   
%   
%   sim = sim.make_excitatory_connection(sim, i, j)
%   
%     Defines a synaptic connection with presynaptic cell 'i' and
%     postsynaptic cell 'j', i.e. 'i' affects 'j'.
%   
%   
%   sim = sim.make_inhibitory_connection(sim, i, j)
%   
%     Defines a synaptic connection with presynaptic cell 'i' and
%     postsynaptic cell 'j', i.e. 'i' affects 'j'.
%   
%   
%   sim = sim.simulate_timeinterval(sim, norecordtime, recordtime)
%   
%     Simulates behavior first for a time given by 'norecordtime' and then
%     for a time given by 'recordtime' during which a record is kept of the
%     membrane potentials for each neuron. 'recordtime' must be positive
%     and nonzero.
%   
%   
%   sim.plot_membranepotentials(sim, figid, neurons)
%   
%     Plots the membrane potentials recorded during a simulation. Function
%     must be called after a simulate_timeinterval() call. 'figid' is the
%     number of the figure window and 'neurons' is a vector of neuron
%     numbers that should be plotted. E.g., plot_membranepotentials(sim,
%     1, [2 4 3]) will plot the records for neurons 2, 4 and 3 in plots
%     appearing in Figure 1.
%   
%   sim = sim.animate_integration(sim, figid, recordtime, frames, neurons)
%   
%     Same as the plotting function above but this one performs the
%     simulation and displays the results 'frames' times

%   January, 2008.
%   Written by Pavel Dimitrov for CPSC475 at Yale University.

sim = make_neuron_simulator(nneurons);

sim.set_Ie = @set_Ie;
sim.simulate_timeinterval = @simulate_timeinterval;
sim.set_mpotential = @set_mpotential;
sim.make_excitatory_connection = @make_excitatory_connection;
sim.make_inhibitory_connection = @make_inhibitory_connection;

sim.plot_membranepotentials = @plot_membranepotentials;
sim.animate_integration = @animate_integration;

%% make_neuron_simulator()
    function sim=make_neuron_simulator(nneurons)
        % function sim=make_neuron_simulator(nneurons)
        sim.nneurons = nneurons;
        sim.Vs = -20e-3*rand(nneurons,1) - 50e-3; % initial membrane potentials
        sim.Ies = 0*ones(nneurons,1); % injected current
        
        sim.PsM = zeros(nneurons, nneurons);     % Matrix of Ps values PsM(i,j) is the Ps in cell j incoming from cell i
        sim.PstimeM = 0.1*ones(nneurons, nneurons); % Matrix of time values for the calculation of Ps(t) (eq. 5.35)
        sim.EsM = zeros(nneurons, nneurons);     % Matrix of synapse current (0mV excitatory and -80mV is inhibitory)
        sim.dt = 1e-4;
    end
%% set_Ie()
    function sim=set_Ie(sim, ni, current)
    % sim -- simulator
    % ni  -- neuron number
    % current -- new injected current value (nA)
    
    sim.Ies(ni) = current*1e-9; 
    end
%% set_mpotential()
    function sim=set_mpotential(sim, ni, mpotential)
    % sim -- simulator
    % ni  -- neuron number
    % mpotential -- new membrane potential in mV
    
    sim.Vs(ni) = mpotential*1e-3; 
    end

%% make_excitatory_connection()
    function sim=make_excitatory_connection(sim, i,j)
        if i<1 || j<1 || i>sim.nneurons || j>sim.nneurons
            error('i and j must be valid neuron numbers, i.e. between 1 and %d', sim.nneurons)
        end

        sim.PsM(i,j)=1;
        sim.EsM(i,j)=0;  % Volts
    end
%% make_inhibitory_connection()
    function sim=make_inhibitory_connection(sim, i,j)
        if i<1 || j<1 || i>sim.nneurons || j>sim.nneurons
            error('i and j must be valid neuron numbers, i.e. between 1 and %d', sim.nneurons)
        end

        sim.PsM(i,j)=1;
        sim.EsM(i,j)=-80*1e-3; % Volts
    end
%% simulate_timestep()
    function sim=simulate_timestep(sim)
        dt = sim.dt;
        
        Pmax = 5;
        tau_s = 10e-3; % time scale for Ps
        
        Vreset = -80*1e-3;    % V
        Vth =    -54*1e-3;    % V
        Vmax =    10*1e-3;    % V

        E_L = -70*1e-3; % (V) leak voltage
        Rm  =  10*1e6;  % ohms, total membrane resistance
        tau =  20*1e-3; % s, tau=cm*rm; % time constant
        
        sim.Vmax = Vmax;
        
        for ni = 1:sim.nneurons
            V =  sim.Vs(ni);
            Ie = sim.Ies(ni);

            if V == Vmax; % is this an action potential
                V = Vreset;
                
                % action potential -> time-since-AP=0
                sim.PstimeM(ni, sim.PsM(ni,:) ~= 0) = 0;
            else
                synapsecontribution=0;
                if find(sim.PsM)
                    for nbn = find(sim.PsM(:,ni)')
                        Ps = sim.PsM(nbn,ni);
                        Es = sim.EsM(nbn,ni);
                        synapsecontribution=synapsecontribution-( 0.05*Ps*(V-Es) );
                    end
                end
                
                dVdt =  1/tau*( E_L - V + Rm*Ie + synapsecontribution); 
                V = V + dt*dVdt;

                if V > Vth
                    V = Vmax;
                end
            end
            
            sim.Vs(ni) = V;
        end
        
        % now update the Ps and the time-since-action-potential PstimeM
        mask = double(sim.PsM ~= 0);
        sim.PstimeM = sim.PstimeM + dt*mask;
        tM = sim.PstimeM;
        sim.PsM = mask.*(Pmax*tM/tau_s.*exp(1 - tM/tau_s));
        
    end

%% simulate_timeinterval()
    function sim=simulate_timeinterval(sim, etint, tint)
        dt = sim.dt;
        
        if tint < 0
           error('recordtime must be nonzeros and positive')
        end
        
        if etint > 0
            timesteps = round(etint/dt);
            for ti=1:timesteps
                sim=simulate_timestep(sim);
            end
        end
        
        timesteps = round(tint/dt);
        t=zeros(1, timesteps);
        VrecordM = zeros(sim.nneurons, timesteps);
        
        for ti=1:timesteps
            % first record the current membrane potentials for each neuron
            for ni=1:sim.nneurons
                VrecordM(ni,ti) = sim.Vs(ni);
            end

            sim=simulate_timestep(sim);
            t(ti+1)=t(ti)+dt;
        end

        sim.VrecordM = VrecordM;
        sim.timev = t(1:end-1); % the time interval vector
    end

%% animate_integration(sim, figid, tint, frames, neuronlist)
    function sim=animate_integration(sim, figid, tint, frames, neuronlist)
        dt = sim.dt;

        timesteps = round(tint/dt);
        t=zeros(1, timesteps);
        VrecordM = zeros(sim.nneurons, timesteps);

        tstp100 = round(timesteps / frames);
        for ti=1:timesteps
            % first record the current membrane potentials for each neuron
            for ni=1:sim.nneurons
                VrecordM(ni,ti) = sim.Vs(ni);
            end

            sim=simulate_timestep(sim);
            t(ti+1)=t(ti)+dt;
            
            % display the current record
            sim.VrecordM = VrecordM;
            sim.timev = t(1:ni); % the time interval vector

            if mod(ti-1, tstp100)==0
                xmin = 0; xmax = tint;
                plot_mps(sim, figid, neuronlist, VrecordM(:,1:ti), t(1:ti), xmin, xmax);
                getframe;
            end
        end

    end

%% plot_membranepotentials()
    function plot_membranepotentials(sim, figid, neuronlist)
        nlsz = length(neuronlist);
        
        figure(figid)
        for dispi=1:nlsz
            subplot(nlsz,1,dispi)
            t = sim.timev;
            v = sim.VrecordM(neuronlist(dispi), :);
            plot(t*1000,v*1000, '-x')
            xlabel('time (ms)')
            ylabel('membrane potential (mV)')
            title(sprintf('Neuron %d (%d spikes)', neuronlist(dispi),...
                sum(v==sim.Vmax)))
        end
    end

    function plot_mps(sim, figid, neuronlist, VrecordM, t, xmin, xmax)
        nlsz = length(neuronlist);
        
        figure(figid)
        for dispi=1:nlsz
            subplot(nlsz,1,dispi)
            %t = sim.timev;
            v = VrecordM(neuronlist(dispi), :);
            plot(t*1000,v*1000, '-x')
            xlabel('time (ms)')
            ylabel('membrane potential (mV)')
            title(sprintf('Neuron %d (%d spikes)', neuronlist(dispi),...
                sum(v==sim.Vmax)))

            ax = axis;
            ax(1)=1000*xmin; ax(2)=1000*xmax; axis(ax);

        end
        
    end
end