# computes the jacobian and the metric for a given model

import jax.numpy as jnp
import numpy as np
from jax import jacfwd, jacrev, jit

from scipy.linalg import expm
from scipy.linalg import eigh
from numpy.linalg import inv
import benchmarks as bm


class FunctionDynamics:
    def __init__(self, model):
        self.model = model

        x = jnp.ones(self.model.dim, dtype='float64')
        if jnp.sum(jnp.abs(self.model.fdyn(0.0, x) - self.model.fdyn(1.0, x))) > 1e-8:
            # https://github.com/google/jax/issues/47
            raise ValueError("Only time-invariant systems supported currently")
        self._cached_f_jac = jit(jacrev(lambda x: self.model.fdyn(0.0, x)))

    def f_jac_at(self, t, x):
        return jnp.array(self._cached_f_jac(x))

    def metric(self, Fmid, ellipsoids):
        if ellipsoids:
            A1inv = Fmid
            A1 = inv(A1inv)
            M1 = inv(A1inv @ A1inv.T)

            W, v = eigh(M1)

            W = abs(W)  # to prevent nan-errors

            semiAxes = 1 / np.sqrt(W)  # needed to compute volume of ellipse

            # Wm = semiAxes * np.eye(self.model.dim)
            # A1inv = (Wm @ v.T).T
            # A1 = inv(A1inv)
            # M1 = A1.T @ A1

        else:
            A1 = np.eye(Fmid.shape[0])
            M1 = np.eye(Fmid.shape[0])
            semiAxes = np.array([1])

        return M1, A1, semiAxes.prod()


def polar2cart_no_rad(phi):
    sin_polar = jnp.sin(phi)
    cart = jnp.append(jnp.cos(phi), jnp.ones(1)) * jnp.append(jnp.ones(1), sin_polar)
    for i in range(1, jnp.size(phi)):
        cart *= jnp.append(jnp.ones(i + 1), sin_polar[:-i])
    return (
        cart  # rad*polar2cart_no_rad(phi) is the true value of the cartesian coordinate
    )

_jac_polar_cached = None


def jacobian_polar_at(phi):
    global _jac_polar_cached
    if _jac_polar_cached is None:
        _jac_polar_cached = jit(jacfwd(polar2cart_no_rad))
    return _jac_polar_cached(phi)


if __name__ == "__main__":
    fdyn = FunctionDynamics(bm.CartpoleCTRNN())

    print(fdyn.f_jac_at(0, fdyn.model.cx))
