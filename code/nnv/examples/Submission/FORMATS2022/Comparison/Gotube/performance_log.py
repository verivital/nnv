import os
import json

registered_args = {}
logged_stats = {}


def log_args(args):
    global registered_args
    for k, v in dict(args).items():
        registered_args[k] = v


def log_stat(stat):
    global logged_stats
    for k, v in stat.items():
        if not k in logged_stats.keys():
            logged_stats[k] = []
        logged_stats[k].append(v)


def close_log(notes):
    global registered_args
    global logged_stats
    final_dict = {"args": registered_args, "stats": logged_stats, "notes": notes}
    os.makedirs("logged", exist_ok=True)
    for i in range(100000):
        fname = f"logged/stat_{i:04d}.json"
        if not os.path.isfile(fname):
            break
    with open(fname, "w") as f:
        json.dump(final_dict, f)


def create_plot_file(files,argsmain):
    fname = files['output_directory']+str(argsmain.benchmark) + "_" + str(argsmain.time_horizon) + "_" + str(argsmain.time_step) + "_" + str(argsmain.batch_size) + "_" + str(argsmain.radius) + "_" + str(argsmain.gamma) + "_" + str(argsmain.mu) +files['output_file']
    return fname


def write_plot_file(fname, mode, t, cx, rad, M1):
    f = open(fname, mode)
    f.write(str(t) + " ")
    f.write(' '.join(map(str, cx.reshape(-1))) + " ")
    f.write(str(rad) + " ")
    f.write(' '.join(map(str, M1.reshape(-1))) + "\n")
    f.close()


if __name__ == "__main__":
    print("registered_args: ", str(registered_args))
    log_args({"hello": "test"})
    print("registered_args: ", str(registered_args))

