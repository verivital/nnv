import numpy as np
import jax.numpy as jnp

import benchmarks as bm
import stochastic_reachtube as reach
import go_tube
import configparser
import time
from performance_log import log_args
from performance_log import close_log
from performance_log import create_plot_file
from performance_log import write_plot_file
from performance_log import log_stat

from scipy.special import gamma

import argparse

# again, this only works on startup!
from jax.config import config
config.update("jax_enable_x64", True) # JAX enforces single precision, need to assign a new configuration for double precision


if __name__ == "__main__":

    start_time = time.time()

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--score", action="store_true")
    parser.add_argument("--benchmark", default="spiralNL")
    # starting_time, time_step and time_horizon for creating reachtubes
    parser.add_argument("--starting_time", default=0.01, type=float)
    parser.add_argument("--time_step", default=0.01, type=float)
    parser.add_argument("--time_horizon", default=10, type=float)
    # batch-size for tensorization
    parser.add_argument("--batch_size", default=10000, type=int)
    # error-probability
    parser.add_argument("--gamma", default=0.01, type=float)
    # mu as maximum over-approximation
    parser.add_argument("--mu", default=1.5, type=float)
    # choose between hyperspheres and ellipsoids to describe the Reachsets
    parser.add_argument("--ellipsoids", action="store_true")
    # initial radius
    parser.add_argument("--radius", default=None, type=float)

    args = parser.parse_args()
    log_args(vars(args))

    config = configparser.ConfigParser()
    config.read("./config.ini")

    files = config["files"]

    rt = reach.StochasticReachtube(
        model=bm.get_model(args.benchmark, args.radius),
        profile=args.profile,
        mu=args.mu,  # mu as maximum over-approximation
        gamma=args.gamma,  # error-probability
        batch=args.batch_size,
        radius=args.radius,
    )  # reachtube

    timeRange = np.arange(args.starting_time, args.time_horizon + 1e-9, args.time_step)
    timeRange_with_start = np.append(np.array([0]), timeRange)

    volume = jnp.array([rt.compute_volume()])

    # propagate center point and compute metric
    (
        cx_timeRange,
        A1_timeRange,
        M1_timeRange,
        semiAxes_prod_timeRange,
    ) = rt.compute_metric_and_center(
        timeRange_with_start,
        args.ellipsoids,
    )

    fname = create_plot_file(files,args)

    write_plot_file(fname, "w", 0, rt.model.cx, rt.model.rad, M1_timeRange[0, :, :])

    total_random_points = None
    total_gradients = None
    total_initial_points = jnp.zeros((0, rt.model.dim))

    # for loop starting at Line 2
    for i, time_py in enumerate(timeRange):
        print(f"Step {i}/{timeRange.shape[0]} at {time_py:0.2f} s")

        rt.time_horizon = time_py
        rt.time_step = args.time_step
        rt.cur_time = rt.time_horizon

        rt.cur_cx = cx_timeRange[i + 1, :]
        rt.A1 = A1_timeRange[i + 1, :, :]

        #  GoTube Algorithm for t = t_j
        #  while loop from Line 5 - Line 12
        (
            rt.cur_rad,
            prob,
            total_initial_points,
            total_random_points,
            total_gradients,
        ) = go_tube.optimize(
            rt, total_initial_points, total_random_points, total_gradients
        )

        write_plot_file(
            fname, "a", time_py, rt.cur_cx, rt.cur_rad, M1_timeRange[i + 1, :, :]
        )

        volume = jnp.append(volume, rt.compute_volume(semiAxes_prod_timeRange[i + 1]))

        if rt.profile:
            # If profiling is enabled, log some statistics about the GD optimization process
            volumes = {
                "volume": float(volume[-1]),
                "average_volume": float(volume.mean()),
            }
            log_stat(volumes)

    if args.score:
        with open("all_prob_scores.csv", "a") as f:
            # CSV with header benchmark, time-horion, prob, runtime, volume
            f.write(f"{args.benchmark},")
            f.write(f"{args.time_horizon:0.4g},")
            f.write(f"{args.radius:0.4g},")
            f.write(f"{1.0-args.gamma:0.4f},")
            f.write(f"{time.time()-start_time:0.2f},")
            f.write(f"{float(volume.mean()):0.5g}")
            f.write("\n")
    if rt.profile:
        final_notes = {
            "total_time": time.time() - start_time,
            "samples": total_random_points.shape[0],
        }
        close_log(final_notes)
