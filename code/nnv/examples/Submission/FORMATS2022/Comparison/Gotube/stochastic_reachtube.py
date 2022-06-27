# optimization problem solved with vanilla gradient descent

import numpy as np
import jax.numpy as jnp
from jax.experimental.ode import odeint
from jax import vmap, partial, jit

from scipy.special import gamma

# own files
import benchmarks as bm
import polar_coordinates as pol
import dynamics


class StochasticReachtube:
    def __init__(
        self,
        model=bm.CartpoleCTRNN(None),
        time_horizon=10.0,  # time_horizon until which the reachtube should be constructed
        profile=False,
        time_step=0.1,  # ReachTube construction
        h_metric=0.05,  # time_step for metric computation
        h_traces=0.01,  # time_step for traces computation
        max_step_metric=0.00125,  # maximum time_step for metric computation
        max_step_optim=0.1,  # maximum time_step for optimization
        samples=100,  # just for plotting: number of random points on the border of the initial ball
        batch=1,  # number of initial points for vectorization
        axis1=0,  # axis to project reachtube to
        axis2=1,
        atol=1e-10,  # absolute tolerance of integration
        rtol=1e-10,  # relative tolerance of integration
        plot_grid=50,
        mu=1.5,
        gamma=0.01,
        radius=False,
    ):

        self.time_step = min(time_step, time_horizon)
        self.profile = profile
        self.h_metric = min(h_metric, time_step)
        self.h_traces = h_traces
        self.max_step_metric = min(max_step_metric, self.h_metric)
        self.max_step_optim = min(max_step_optim, self.time_step)
        self.time_horizon = time_horizon
        self.samples = samples
        self.batch = batch
        self.axis1 = axis1
        self.axis2 = axis2
        self.atol = atol
        self.rtol = rtol
        self.plotGrid = plot_grid
        self.mu = mu
        self.gamma = gamma

        self.model = model
        self.init_model()

        self.metric = dynamics.FunctionDynamics(model).metric
        self.init_metric()

        self.f_jac_at = dynamics.FunctionDynamics(model).f_jac_at

    def init_metric(self):
        self.M1 = np.eye(self.model.dim)
        self.A1 = np.eye(self.model.dim)
        self.A1inv = np.eye(self.model.dim)
        self.A0inv = np.eye(self.model.dim)

    def init_model(self):
        self.cur_time = 0
        self.cur_cx = self.model.cx
        self.cur_rad = self.model.rad
        self.t0 = 0
        self.cx_t0 = self.model.cx
        self.rad_t0 = self.model.rad

    def compute_volume(self, semiAxes_product=None):
        if semiAxes_product is None:
            semiAxes_product = 1
        volC = gamma(self.model.dim / 2.0 + 1) ** -1 * jnp.pi ** (
            self.model.dim / 2.0
        )  # volume constant for ellipse and ball
        return volC * self.cur_rad ** self.model.dim * semiAxes_product

    def plot_traces(self, axis_3d):
        rd_polar = pol.init_random_phi(self.model.dim, self.samples)
        rd_x = (
            vmap(pol.polar2cart, in_axes=(None, 0))(self.model.rad, rd_polar)
            + self.model.cx
        )
        plot_timerange = jnp.arange(0, self.time_horizon + 1e-9, self.h_traces)

        sol = odeint(
            self.fdyn_jax,
            rd_x,
            plot_timerange,
            atol=self.atol,
            rtol=self.rtol,
        )

        for s in range(self.samples):
            axis_3d.plot(
                xs=sol[:, s, self.axis1],
                ys=sol[:, s, self.axis2],
                zs=plot_timerange,
                color="k",
                linewidth=1,
            )

        p_dict = {
            "xs": np.array(sol[:, s, self.axis1]),
            "ys": np.array(sol[:, s, self.axis2]),
            "zs": np.array(plot_timerange),
        }
        return p_dict

    def propagate_center_point(self, time_range):
        cx_jax = self.model.cx.reshape(1, self.model.dim)
        F = jnp.eye(self.model.dim)
        aug_state = jnp.concatenate((cx_jax, F)).reshape(1, -1)
        sol = odeint(
            self.aug_fdyn_jax,
            aug_state,
            time_range,
            atol=self.atol,
            rtol=self.rtol,
        )
        cx, F = vmap(self.reshape_aug_state_to_matrix)(sol)
        return cx, F

    def compute_metric_and_center(self, time_range, ellipsoids):
        cx_timeRange, F_timeRange = self.propagate_center_point(time_range)
        A1_timeRange = np.eye(self.model.dim).reshape(1, self.model.dim, self.model.dim)
        M1_timeRange = np.eye(self.model.dim).reshape(1, self.model.dim, self.model.dim)
        semiAxes_prod_timeRange = np.array([1])

        for idx, t in enumerate(time_range[1:]):
            M1_t, A1_t, semiAxes_prod_t = self.metric(
                F_timeRange[idx + 1, :, :], ellipsoids
            )
            A1_timeRange = np.concatenate(
                (A1_timeRange, A1_t.reshape(1, self.model.dim, self.model.dim)), axis=0
            )
            M1_timeRange = np.concatenate(
                (M1_timeRange, M1_t.reshape(1, self.model.dim, self.model.dim)), axis=0
            )
            semiAxes_prod_timeRange = np.append(
                semiAxes_prod_timeRange, semiAxes_prod_t
            )

        return cx_timeRange, A1_timeRange, M1_timeRange, semiAxes_prod_timeRange

    def reshape_aug_state_to_matrix(self, aug_state):
        aug_state = aug_state.reshape(-1, self.model.dim)  # reshape to matrix
        x = aug_state[:1][0]
        F = aug_state[1:]
        return x, F

    def reshape_aug_fdyn_return_to_vector(self, fdyn_return, F_return):
        return jnp.concatenate((jnp.array([fdyn_return]), F_return)).reshape(-1)

    @partial(jit, static_argnums=(0,))
    def aug_fdyn(self, t=0, aug_state=0):
        x, F = self.reshape_aug_state_to_matrix(aug_state)
        fdyn_return = self.model.fdyn(t, x)
        F_return = jnp.matmul(self.f_jac_at(t, x), F)
        return self.reshape_aug_fdyn_return_to_vector(fdyn_return, F_return)

    def aug_fdyn_jax(self, aug_state=0, t=0):
        return vmap(self.aug_fdyn, in_axes=(None, 0))(t, aug_state)

    def fdyn_jax(self, x=0, t=0):
        return vmap(self.model.fdyn, in_axes=(None, 0))(t, x)

    def create_aug_state(self, polar, rad_t0, cx_t0):
        x = jnp.array(
            pol.polar2cart_euclidean_metric(rad_t0, polar, self.A0inv) + cx_t0
        )
        F = jnp.eye(self.model.dim)

        aug_state = jnp.concatenate((jnp.array([x]), F)).reshape(
            -1
        )  # reshape to row vector

        return aug_state, x

    def create_aug_state_cartesian(self, x, F):
        aug_state = jnp.concatenate((jnp.array([x]), F)).reshape(
            -1
        )  # reshape to row vector

        return aug_state

    def one_step_aug_integrator(self, x, F):
        aug_state = vmap(self.create_aug_state_cartesian, in_axes=(0, 0))(x, F)
        sol = odeint(
            self.aug_fdyn_jax,
            aug_state,
            jnp.array([0, self.time_step]),
            atol=self.atol,
            rtol=self.rtol,
        )
        x, F = vmap(self.reshape_aug_state_to_matrix)(sol[-1])
        return x, F

    def aug_integrator(self, polar, step=None):
        if step is None:
            step = self.cur_time

        rad_t0 = self.rad_t0
        cx_t0 = self.cx_t0

        aug_state, initial_x = vmap(self.create_aug_state, in_axes=(0, None, None))(
            polar, rad_t0, cx_t0
        )
        sol = odeint(
            self.aug_fdyn_jax,
            aug_state,
            jnp.array([0, step]),
            atol=self.atol,
            rtol=self.rtol,
        )
        x, F = vmap(self.reshape_aug_state_to_matrix)(sol[-1])
        return x, F, initial_x

    def aug_integrator_neg_dist(self, polar):
        x, F, initial_x = self.aug_integrator(polar)
        neg_dist = vmap(self.neg_dist_x)(x)
        return x, F, neg_dist, initial_x

    def one_step_aug_integrator_dist(self, x, F):
        x, F = self.one_step_aug_integrator(x, F)
        neg_dist = vmap(self.neg_dist_x)(x)
        return x, F, -neg_dist

    def neg_dist_x(self, xt):
        dist = jnp.linalg.norm(jnp.matmul(self.A1, xt - self.cur_cx))
        return -dist