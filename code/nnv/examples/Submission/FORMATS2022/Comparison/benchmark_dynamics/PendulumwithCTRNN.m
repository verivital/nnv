function dx = PendulumwithCTRNN(x,u)
% class PendulumwithCTRNN:
%     def __init__(self, radius):
%         # ============ adapt initial values ===========
%         if radius is not None:
%             self.rad = radius
%         else:
%             self.rad = 0.5
%         # ===================================================
%         self.cx = np.zeros(10)
%         self.dim = self.cx.size  # dimension of the system
%         arr = np.load("rl/pendulum_ctrnn.npz")
%         self.params = {k: arr[k] for k in arr.files}
% 
%     def fdyn(self, t=0, x=None):
%         if x is None:
%             x = np.zeros(self.dim, dtype=object)

        w1 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/w1.npy');
        w2 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/w2.npy');
        wa = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/wa.npy');
        b1 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/b1.npy');
        b2 = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/b2.npy');
        ba = readNPY('/home/manzand/Documents/MATLAB/neuralODE/Gotube/benchmark_dynamics/rl/lds/ba.npy');

%         hidden = np.tanh(np.dot(x, self.params["w1"]) + self.params["b1"])
        hidden = tanh(w1' * x + b1);
%         dhdt = 'np.dot(hidden, self.params["w2"]) + self.params["b2"]
        dhdt = w2' * x + b2;
% 
%         action = np.tanh(np.dot(hidden, self.params["wa"]) + self.params["ba"])
        action = tanh(wa' * hidden + ba);
%         th, thdot = x[-2:]
        th = x(9);
% 
        max_speed = 8;
        g = 9.81;
        m = 1.0;
        l = 1.0;
% 
        newthdot = -3 * g / (2 * l) * sin(th + pi) + 3.0 / (m * l ^ 2) * action;
        newthdot = max_speed * tanh(newthdot / max_speed);
        newth = newthdot;
% 
%         dxdt = np.array([newth]).reshape((1,))
        dxdt = newth;
%         dydt = np.array([newthdot]).reshape((1,))
        dydt = newthdot;
%         dfdt = np.concatenate([dhdt, dxdt, dydt], axis=0)
%         return dfdt
        dx = [dhdt; dxdt; dydt];

end

