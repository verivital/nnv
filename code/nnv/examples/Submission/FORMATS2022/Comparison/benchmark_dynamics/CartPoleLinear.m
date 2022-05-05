function dx = CartPoleLinear(x,u)
% # 4-dimensional cartpole with linear stabilizing controller
% class CartpoleLinear:
%     def __init__(self, radius):
%         # ============ adapt initial values ===========
%         self.cx = (0, 0, 0.001, 0)  # initial values
%         if radius is not None:
%             self.rad = radius
%         else:
%             self.rad = 0.05
%         # ===================================================
% 
%         self.cx = np.array(self.cx, dtype=float)
%         self.dim = self.cx.size  # dimension of the system
% 
%     def fdyn(self, t=0, x=None):
%         if x is None:
%             x = np.zeros(self.dim, dtype=object)
% 
%         # ============ adapt input and system dynamics ===========
%         dth, dx, th, x = x  # input variables
% 
%         M = 1.0
%         g = 9.81
%         l = 1.0
%         m = 0.001
% 
%         f = -1.1 * M * g * th - dth
% 
%         fdth = (
%             1.0
%             / (l * (M + m * sin(th) * sin(th)))
%             * (
%                 f * cos(th)
%                 - m * l * dth * dth * cos(th) * sin(th)
%                 + (m + M) * g * sin(th)
%             )
%         )
%         fdx = (
%             1.0
%             / (M + m * sin(th) * sin(th))
%             * (f + m * sin(th) * (-l * dth * dth + g * cos(th)))
%         )
% 
%         fx = dx
% 
%         fth = dth
% 
%         system_dynamics = [
%             fdth,
%             fdx,
%             fth,
%             fx,
%         ]  # has to be in the same order as the input variables
%         # ===================================================
% 
%         return np.array(system_dynamics)  # return as numpy array

% x1 = dth; x2 = dx; x3 = th; x4 = x;

M = 1.0;
g = 9.81;
l = 1.0;
m = 0.001;

f = -1.1 * M * g * x(3) - x(1);

dx(1,1) = (1.0 / (l * (M + m * sin(x(3)) * sin(x(3)))) * (f * cos(x(3)) - m * l * x(1) * x(1) * cos(x(3)) * sin(x(3))+ (m + M) * g * sin(x(3))));

dx(2,1) = (1.0 / (M + m * sin(x(3)) * sin(x(3))) * (f + m * sin(x(3)) * (-l * x(1) * x(1) + g * cos(x(3)))));

dx(3,1) = x(1);

dx(4,1) = x(4);

end

