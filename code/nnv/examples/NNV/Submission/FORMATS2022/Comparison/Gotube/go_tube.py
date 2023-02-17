# Algorithms of SLR paper for safety region, probability and stoch. optimization

import jax.numpy as jnp
from jax import vmap
import polar_coordinates as pol
from jax.numpy.linalg import svd
import jax.scipy.special as sc
import time
from performance_log import log_stat
from timer import Timer
from scipy import stats


# using expected difference quotient of center lipschitz constant
def get_safety_region_radius(model, dist, dist_best, lip, lip_mean_diff):
    safety_radius = -lip + jnp.sqrt(lip ** 2 + 4 * lip_mean_diff * (model.mu * dist_best - dist))
    safety_radius = safety_radius / (2 * lip_mean_diff)
    safety_radius = safety_radius * (dist > 0) * (lip_mean_diff > 0)

    safety_radius = jnp.minimum(safety_radius, 2 * model.rad_t0)

    return safety_radius


def compute_maximum_singular_value(model, F):
    F_metric = jnp.matmul(model.A1, F)
    F_metric = jnp.matmul(F_metric, model.A0inv)
    _, sf, _ = svd(F_metric)
    max_sf = jnp.max(sf)

    return max_sf


def get_angle_of_cap(model, radius):
    radius = jnp.minimum(radius, 2 * model.rad_t0)
    return 2 * jnp.arcsin(0.5 * radius / model.rad_t0)


def get_probability_of_cap(model, radius):
    with Timer('get angle of cap'):
        angle = get_angle_of_cap(model, radius)
    with Timer('get probability of cap'):
        a = 0.5 * (model.model.dim - 1)
        b = 0.5
        x = jnp.sin(angle) ** 2
        betainc_angle = 0.5 * sc.betainc(a, b, x)

        # formula is only for the smaller cap with angle <= pi/2, sinus is symmetric => thus use 1-area otherwise
        betainc_angle = jnp.where(angle > 0.5 * jnp.pi, 1 - betainc_angle, betainc_angle)

    return betainc_angle


def get_probability_not_in_cap(model, radius):
    return 1 - get_probability_of_cap(model, radius)


def get_probability_none_in_cap(model, radius_points):
    return jnp.prod(get_probability_not_in_cap(model, radius_points))


#  probability calculation using http://docsdrive.com/pdfs/ansinet/ajms/2011/66-70.pdf (equation 1
#  page 68) and the normalized incomplete Beta-Function in scipy (
#  https://scipy.github.io/devdocs/generated/scipy.special.betainc.html#scipy.special.betainc) - Only use the
#  random sampled points for probability construction
#  use also the discarded points and create balls around them
def get_probability(model, radius_points):
    return jnp.sqrt(1-model.gamma) * (1 - get_probability_none_in_cap(model, radius_points))


def get_diff_quotient(x, fx, y_jax, fy_jax, axis):
    distance = jnp.linalg.norm(x - y_jax, axis=axis)
    diff_quotients = abs(fx - fy_jax) / distance * (distance > 0)
    return diff_quotients


def get_diff_quotient_vmap(x_jax, fx_jax, y_jax, fy_jax, axis):
    return vmap(get_diff_quotient, in_axes=(0, 0, None, None, None))(x_jax, fx_jax, y_jax, fy_jax, axis)


def optimize(model, initial_points, previous_points=None, previous_gradients=None):
    start_time = time.time()

    prob = None

    sample_size = model.batch
    df = sample_size - 2
    conf = (1 + jnp.sqrt(1-model.gamma)) / 2
    t_star = stats.t.ppf(conf, df)

    if previous_points is None or previous_gradients is None:
        previous_samples = 0
        phis = pol.init_random_phi(model.model.dim, model.batch)
        points, gradients, neg_dists, initial_points = model.aug_integrator_neg_dist(phis)
        dists = -neg_dists
    else:
        previous_samples = previous_points.shape[0]
        with Timer('integrate random points and gradients - one step'):
            points, gradients, dists = model.one_step_aug_integrator_dist(
                previous_points, previous_gradients
            )

    first_iteration = True

    while prob is None or prob < 1 - model.gamma:

        if not first_iteration:
            with Timer('sample phis'):
                phis = pol.init_random_phi(model.model.dim, model.batch)
            with Timer('compute first integration step and dist'):
                new_points, new_gradients, new_neg_dists, new_initial_points = model.aug_integrator_neg_dist(phis)
                new_dists = -new_neg_dists

            with Timer('concatenate new points to tensors'):
                points = jnp.concatenate((points, new_points), axis=0)
                gradients = jnp.concatenate((gradients, new_gradients), axis=0)
                dists = jnp.concatenate((dists, new_dists), axis=0)
                initial_points = jnp.concatenate((initial_points, new_initial_points), axis=0)

        with Timer('compute best dist'):
            dist_best = dists.max()

        with Timer('compute lipschitz'):
            # compute maximum singular values of all new gradient matrices
            lipschitz = vmap(compute_maximum_singular_value, in_axes=(None, 0))(model, gradients)

        with Timer('compute expected local lipschitz'):
            # compute expected value of delta lipschitz
            dimension_axis = 1
            # limit expected value to batch size
            diff_quotients = get_diff_quotient_vmap(
                initial_points,
                lipschitz,
                initial_points[:sample_size],
                lipschitz[:sample_size],
                dimension_axis
            )

            v_mean = jnp.nanmean(diff_quotients, axis=dimension_axis)
            v_std = jnp.nanstd(diff_quotients, axis=dimension_axis)

            delta_lipschitz = v_mean + t_star * v_std / jnp.sqrt(sample_size)

        with Timer('get safety region radii'):
            safety_region_radii = get_safety_region_radius(
                model, dists, dist_best, lipschitz, delta_lipschitz
            )

        with Timer('compute probability'):
            prob = get_probability(model, safety_region_radii)

        if first_iteration:
            print("start probability is: ")
            print(prob)
        else:
            print("current probability is: ")
            print(prob)
            print("number of samples: ")
            print(points.shape[0])

        first_iteration = False

    print('prob after loop: %s' % prob)

    new_samples = points.shape[0] - previous_samples

    print(
        f"Visited {new_samples} new points in {time.time() - start_time:0.2f} seconds."
        # Current probability coverage {100.0 * prob:0.3f}%"
    )
    print("Probability reached given value!")

    dist_with_safety_mu = model.mu * dist_best

    if model.profile:
        # If profiling is enabled, log some statistics about the GD optimization process
        stat_dict = {
            "loop_time": time.time() - start_time,
            "new_points": int(new_samples),
            "total_points": int(previous_samples + new_samples),
            "prob": float(prob),
            "dist_best": float(dist_best),
            "radius": float(dist_with_safety_mu),
        }
        log_stat(stat_dict)

    return dist_with_safety_mu, prob, initial_points, points, gradients
