clear;
load controller_3_20.mat;
weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{i, 1}, 'ReLU');
    Layers = [Layers L];
end

L = Layer(weights{1, n}, bias{n, 1}, 'Linear');

Layers = [Layers L];

Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(6, 1, @car_dynamics);
Plant.set_timeStep(0.01); % time step for reachability analysis of the plant
Plant.set_tFinal(0.1); % Ts = 0.1, sampling time for control signal from neural network controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant

% x = [x_lead v_lead x_internal_lead x_ego v_ego x_internal_ego]'

lb = [65; 32; 0; 10; 30; 0];
ub = [70; 32.2; 0; 11; 30.2; 0];

init_set = Box(lb, ub);

% reference input for neural network controller
t_gap = 1.4;
D_default = 10;

lb1 = [30; 1.4];
ub1 = [30; 1.4];

input_ref = Box(lb1, ub1);

N = 50;  

n_cores = 4; % number of cores 

% safety specification dis > safe_dis

unsafe_mat = [1 0 0 -1 -t_gap 0];
unsafe_vec = [D_default];

n_samples = 10; 

[falsify_result, falsify_time, counter_sim_traces, counter_control_traces, counter_init_states, counter_ref_inputs] = ncs.falsify(0.1, N, init_set, input_ref, unsafe_mat, unsafe_vec, n_samples);

if falsify_result == 1
    
    for i=1:length(counter_sim_traces)
        simTrace = cell2mat(counter_sim_traces{1, i});
        controlTrace = cell2mat(counter_control_traces{1, i});
        dis = [1 0 0 -1 0 0] * simTrace; 
        safe_dis = [0 0 0 0 t_gap 0] * simTrace + D_default;
        times = 0:0.1:0.1*N;
        figure;
        % subplot 1
        subplot(3, 1, 1);

        plot(times, simTrace(2, :), '-', 'Color', 'green');
        hold on;
        plot(times, simTrace(5, :), '-', 'Color', 'blue');
        xlabel('time (seconds)');
        ylabel('Velocity');
        title('Ego car velocity (blue) vs. lead car velocity (green)');
        % subplot 2
        subplot(3,1,2);
        plot(times, dis, '-', 'Color', 'blue');
        hold on;
        plot(times, safe_dis, '-', 'Color', 'red');
        xlabel('time (seconds)');
        ylabel('Distance');
        title('Relative distance (blue) vs. safe distance (red)');
        % subplot 3
        subplot(3,1,3);
        plot(times(1:N), controlTrace, '-', 'Color', 'cyan');
        xlabel('time (seconds)');
        ylabel('Control')
        title('Control vs. time');
    end
    

    saveas(gcf, 'falsifyTrace.pdf');
end

