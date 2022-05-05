function dx = LDSwithCTRNN(x,u)
% class LDSwithCTRNN:
%     def __init__(self, radius):
%         # ============ adapt initial values ===========
%         if radius is not None:
%             self.rad = radius
%         else:
%             self.rad = 0.5
%         # ===================================================
%         self.cx = np.zeros(10)
%         self.dim = self.cx.size  # dimension of the system
%         arr = np.load("rl/lds_ctrnn.npz")
%         self.params = {k: arr[k] for k in arr.files}
% 
%     def fdyn(self, t=0, x=None):
%         if x is None:
%             x = np.zeros(self.dim, dtype=object)
        % Load parameters
        w1 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/w1.npy');
        w2 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/w2.npy');
        wa = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/wa.npy');
        b1 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/b1.npy');
        b2 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/b2.npy');
        ba = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/ba.npy');

        hidden = tanh(w1' * x + b1); % tanh
        dhdt = w2' * hidden + b2; % linear
% 
        action = tanh(wa' * hidden + ba); % tanh
%         x, y = x[-2:]
%         x2 = x(2);
%         y = x(3);
% 
        dxdt = x(10);
        dydt = 0.2 + 0.4 * action;
% 
%         dxdt = np.array([dxdt]).reshape((1,))
%         dydt = np.array([dydt]).reshape((1,))
%         dfdt = np.concatenate([dhdt, dxdt, dydt], axis=0)
        dx = [dhdt;dxdt;dydt];
end

