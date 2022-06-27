# plotting the outputs from LRT-NG (ellipse, circle, intersections etc.)

import numpy as np
import matplotlib.pyplot as plt

from numpy.linalg import svd
from numpy.linalg import inv

import stochastic_reachtube as reach
import benchmarks as bm
import configparser
import argparse
import pickle
import os


def draw_ellipse(ellipse, color, alph, axis_3d):
    c = ellipse[1:3]
    r = ellipse[3]
    M = np.reshape(ellipse[4:8], (2, 2)) / r ** 2

    # "singular value decomposition" to extract the orientation and the
    # axes of the ellipsoid
    _, D, V = svd(M)

    plot_grid = 50

    # get the major and minor axes
    a = 1 / np.sqrt(D[0])
    b = 1 / np.sqrt(D[1])

    theta = np.arange(0, 2 * np.pi + 1 / plot_grid, 1 / plot_grid)

    # parametric equation of the ellipse
    state = np.zeros((2, np.size(theta)))
    state[0, :] = a * np.cos(theta)
    state[1, :] = b * np.sin(theta)

    # coordinate transform
    X = V @ state

    X[0] += c[0]
    X[1] += c[1]

    axis_3d.plot(xs=X[0], ys=X[1], zs=ellipse[0], color=color, alpha=alph)


def plot_ellipse(
    time_horizon, dim, axis1, axis2, file, color, alph, axis_3d, skip_reachsets=1
):
    data_ellipse = np.loadtxt(file)

    # permutation matrix to project on axis1 and axis2
    P = np.eye(dim)
    P[:, [0, axis1]] = P[:, [axis1, 0]]
    P[:, [1, axis2]] = P[:, [axis2, 1]]

    count = skip_reachsets

    for ellipse in data_ellipse[1:]:

        if count != skip_reachsets:
            count += 1
            continue

        count = 1

        ellipse2 = ellipse

        # create ellipse plotting values for 2d projection
        # https://math.stackexchange.com/questions/2438495/showing-positive-definiteness-in-the-projection-of-ellipsoid
        # construct ellipse2 to have a 2-dimensional ellipse as an input to
        # ellipse_plot

        if dim > 2:
            center = ellipse[1 : dim + 1]
            ellipse2[1] = center[axis1]
            ellipse2[2] = center[axis2]
            radius_ellipse = ellipse[dim + 1]
            m1 = np.reshape(ellipse[dim + 2 :], (dim, dim))
            m1 = m1 / radius_ellipse ** 2
            m1 = P.transpose() @ m1 @ P  # permutation to project on chosen axes
            ellipse2[3] = 1  # because radius is already in m1

            # plot ellipse onto axis1-axis2 plane
            J = m1[0:2, 0:2]
            K = m1[2:, 2:]
            L = m1[2:, 0:2]
            m2 = J - L.transpose() @ inv(K) @ L
            ellipse2[4:8] = m2.reshape(1, -1)

        draw_ellipse(ellipse2[0:8], color, alph, axis_3d)

        if ellipse[0] >= time_horizon:
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--time_step", default=0.01, type=float)
    parser.add_argument("--time_horizon", default=0.01, type=float)
    parser.add_argument("--benchmark", default="bruss")
    parser.add_argument("--output_number", default="0000")
    parser.add_argument("--samples", default=100, type=int)
    parser.add_argument("--axis1", default=0, type=int)
    parser.add_argument("--axis2", default=1, type=int)
    # initial radius
    parser.add_argument("--radius", default=None, type=float)

    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read("config.ini")

    files = config["files"]

    rt = reach.StochasticReachtube(
        model=bm.get_model(args.benchmark, args.radius),
        time_horizon=args.time_horizon,
        time_step=args.time_step,
        samples=args.samples,
        axis1=args.axis1,
        axis2=args.axis2,
    )  # reachtube

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    p_dict = rt.plot_traces(axis_3d=ax)
    p_dict["time_horizon"] = args.time_horizon
    p_dict["dim"] = rt.model.dim
    p_dict["axis1"] = rt.axis1
    p_dict["axis2"] = rt.axis2
    p_dict["data_ellipse"] = np.loadtxt(
        files["output_directory"] + str(args.output_number) + files["output_file"]
    )

    plot_ellipse(
        args.time_horizon,
        rt.model.dim,
        rt.axis1,
        rt.axis2,
        files["output_directory"] + str(args.output_number) + files["output_file"],
        color="magenta",
        alph=0.8,
        axis_3d=ax,
        skip_reachsets=1,
    )

    plt.show()

    os.makedirs("plot_obj", exist_ok=True)
    for i in range(1000):
        filename = f"plot_obj/plot_{i:03d}.pkl"
        if not os.path.isfile(filename):
            with open(filename, "wb") as f:
                pickle.dump(p_dict, f, protocol=pickle.DEFAULT_PROTOCOL)
            break
