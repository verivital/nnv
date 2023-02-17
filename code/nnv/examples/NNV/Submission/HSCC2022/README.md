This code allows verifying safety properties of a particular class of system, based upon the *Neural Network Verification toolbox* (NNV) and the Github repository *Verification and simulation of AcasXu closed-loop system* [1].
The type of system that can be analyzed has the following caracteristics:
- it is composed of multiple (possibly one) cyber-physical systems e.g., a fleet of UAVs
- each CPS is equipped with a controller that is a classifier based on multiple neural networks

This type of system is further described in the paper *Verification of machine learning based cyber-physical systems: a comparative study*, accepted to HSCC 2022.

The code provides:
1. the models for three systems:
	- the Vertical Collision Avoidance System (VCAS) which is the combination of two aircraft, each one equipped with a controller providing vertical collision avoidance advisories. The associated safety criteria is the absence of a collision between the two aircraft.
	- the Airborne Collision Avoidance System (ACAS Xu) which is the combination of two aircraft and which has two variants: ACASXu2Agent where the two aircraft are equipped with a controller providing horizontal avoidance advisories and ACASXu1Agent where one of the two aircraft follows a uniform rectilinear trajectory. The associated safety criteria is the same as for the VCAS.
	- the Cartpole, which is an inverted pendulum placed on a moving cart, equipped with a controller that prescribes the displacement of the cart. The system is safe when the cartpole remains balanced i.e., the inverted pendulum does not fall.

    A detailed description of these three systems is given in the paper cited above.

2. an API to run any user-defined verification problem. A verification problem is defined as:
	- a system (it can be one of the three systems mentioned above but one could define another system: the code organization is described below)
	- the initial state of the physical component of the system
	- an uncertainty on this initial state
	- the initial state of the controller(s) i.e., the initial command(s), produced at the previous control step
	- a time horizon

    See below for instructions on how to run such a problem.


## Installation

Follow the installation guidelines of NNV.


## Run a user-defined verification problem

Open Matlab, go to the `nnv/code/nnv` directory and run the `install.m` and `startup_nnv.m` scripts.

Then go tho this directory and type the following commands in the console of Matlab:

1. Load the system

    `>> Param = SetUpVCAS2AgentVerif; % example for the VCAS`

    NOTE: 
    - use `SetUpVCAS2AgentVerif` for the VCAS system
    - use `SetUpACASXu1AgentVerif` for the ACASXu1Agent system
    - use `SetUpACASXu2AgentVerif` for the ACASXu2Agent system
    - use `SetUpCartPoleVerif for` the Cartpole system

2. Specify the value of the time horizon t_end

    ``>> t_end = 5.0;``

3. Specify the initial state of the pysical component of the system in the format `[var_state_1,var_state_2,...,var_state_m]`

    ``>> init_state = [-100.0,5.0,-5.0,5.0]; % example for the VCAS system``

    NOTE:
    - for the VCAS, the initial state is in the following format: $[h, \dot{h_\text{own}},\dot{h_\text{int}},\tau]$
    - for the ACAS, the initial state is in the following format: $[x_\text{own},y_\text{own},\pi/2+\psi_\text{own},x_\text{int},y_\text{int},\pi/2+\psi_\text{int},v_\text{own},v_\text{int}]$
    - for the Cartole, the initial state is in the following format: $[x,\dot{x},\theta,\dot{\theta}]$

4. Specify the initial command(s) in the format `[idx_init_command_CPS_1,...,idx_init_command_CPS_n]`

    ``>> init_command = [1,1]; % example for the VCAS system``

    NOTE:
    - for the VCAS, the initial command is in the format `[idx_1,idx_2]` where both 
idx_1 and idx_2 range from 1 to 9 (9 possible commands for each aircraft)
    - for the ACASXu2Agent, the initial command is in the format `[idx_1,idx_2]` where both 
idx_1 and idx_2 range from 1 to 5 (5 possible commands for each aircraft)
    - for the ACASXu1Agent, the initial command is in the format `[idx_1]` where idx_1 
ranges from 1 to 5 (5 possible commands for the ownship)
    - for the Cartpole, the initial command is in the format `[idx_1]` where idx_1 
ranges from 1 to 2 (2 possible commands for the cartpole)

5. Specify the uncertainty on the initial state

    ``>> uncertainty = [5.0,0.0,0.0,0.0]; % example for the VCAS system``

6. Build the initial set, based on the initial state and the uncertainty

    ``>> init_set = InitializeState(Param, init_state, uncertainty);``

7. Run the safety verification program

    `>> [allReach, verifRes, t_elapsed, t_dynamics_eval, t_nn_eval, max_n_branches] = ClosedLoopVerif(Param, init_set, init_command, t_end);`

    Then you can display the verification result (safe or unknown) by typing:
    
    ``>> verifRes``

    You can also display the verification time by typing:
    
    ``>> t_elapsed``
    
    The allReach cell struct contains the reachable sets over time.
    
8. Plot the results of the reachability analysis

    ``>> Param.plotResults(allReach, init_state);``
    

## Code organization


The code contains two types of functions:

1. Generic functions:

	- `InitializeState` builds the set of the possible initial states of the physical component of the system, by considering a given uncerainty on a given initial state.
	- `ClosedLoopVerif` constructs a flowpipe that over-approximates the reachable states of the system, starting from a set of possible initial states. It alternates between approximating the dynamics of the physical component of the system and the behaviour of the controller. To this end, it calls:
        - The `PlantReach` function that uses the CORA tool (with zonotope based approximation) to approximate the dynamics of the physical component of the system
	    - The functions related to the controller, among which:
            - The `Normalize` function that normalizes the intputs of the neural networks
            - The `ReachNN` function that analyzes a neural network with approximated star set abstraction
            - The `PostProcessing` function that computes the reachable commands associated to the output set of a neural network (this function uses the `ArgMaxVerif` or `ArgMinVerif` function)

    The code also contains some auxiliary functions such as `CreateNNStruct`, `LimitAngleSet_PiPi`, `preSetCombos`, `subCombos`.

2. System specific functions: the system specific features are specified through the `SetUp` functions (*e.g.,* `SetUpVCAS2AgentVerif` for the VCAS).
This function points to additional materials (code and neural networks) which are located in the `Models` directory.

    For example, the VCAS specific code is located in the `Models/VCAS` directory. It contains:

    - The ReLU neural networks used by the controller, stored as `.mat` files in the `nnv_format` directory
    - The `VCASVerifNeuralNetworks` function that loads the neural networks
    - The `VCASVerifPreProcessing` function, that represents the pre-processing function of the controller
    - The `VCASVerifModel2Agent` function that represents the differential equation driving the dynamics of the physical system
    - The `VCASVerifSafetyness` function that checks whether a given state of the CPS is safe or not
    - The `VCASVerifStoppingCriteria` function which can stop the analysis if a given stopping criteria is reached, for example the system has terminated its mission
    - The `VCASPlot` function which plots the results of the reachability analysis


    
[1] https://github.com/mldiego/AcasXu
