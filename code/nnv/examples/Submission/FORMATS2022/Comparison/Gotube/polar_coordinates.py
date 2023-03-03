# transformation between polar and cartesian coordinates

import numpy as np
import jax.numpy as jnp
from jax import jit, partial
import dynamics

# initialize random polar coordinates with dimension dim
# TODO @Mathias please check if this is fine like that for random process
_rng = np.random.RandomState(12937)


@partial(jit, static_argnums=(0, 1))
def init_random_phi(dim, samples=1):
    global _rng
    phi = _rng.uniform(0, jnp.pi, samples * (dim - 2))
    phi = jnp.append(phi, _rng.uniform(0, 2 * jnp.pi, samples))
    phi = jnp.reshape(phi, (samples, dim - 1), order="F")

    return phi


@jit
def polar2cart(rad, phi):
    return rad * polar2cart_no_rad(phi)


def polar2cart_euclidean_metric(rad, phis, A0inv):
    return rad * jnp.matmul(A0inv, polar2cart_no_rad(phis))


def polar2cart_no_rad(phi):
    return dynamics.polar2cart_no_rad(phi)