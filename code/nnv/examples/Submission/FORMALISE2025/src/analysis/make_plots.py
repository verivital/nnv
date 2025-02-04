import csv
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

# local modules
import vvn.prep as vp
from vvn.config import Config
import vvn.gtsrbprep as vgp
import vvn.stmnistprep as vsp

x = ['1/255', '2/255', '3/255']

# get the results dir
root = os.path.dirname(os.getcwd()) # Submission folder
output_dir = os.path.join(root, 'FORMALISE2025' ,'results')

def summarize(output_file_dir):
    average_times = []
    for filename in os.listdir(output_file_dir):
        if filename == '.DS_Store':
            continue

        fp = os.path.join(output_file_dir, filename)

        # open the results csv file
        with open(fp, 'r', newline='') as f:
            reader = csv.reader(f, delimiter=',')

            # skip the header
            next(reader)

            res = []
            t = []

            # read the values and build the new arrays for analysis
            for row in reader:
                res.append(row[1])
                t.append(row[2] if not row[2] == 'timeout' else 1800.0)

            # have to convert strings to valid floats before casting to int
            res = np.array(res).astype(float)
            res = np.array(res).astype(int)
            t = np.array(t).astype(float)

        # count the number of verified samples
        total_verified = np.sum(res[res == 1])

        # calculate average time to verify
        average_time = np.mean(t)

        average_times.append(average_time)

    return average_times

# ZOOM IN
config = Config(
        class_size=10,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='zoom_in',
        sample_len=4,
        ver_algorithm='relax',
        timeout=1800,
        output_dir=output_dir
)
zoomin_4 = summarize(vp.build_output_filepath(config=config, parent_only=True))

config.sample_len = 8
zoomin_8 = summarize(vp.build_output_filepath(config=config, parent_only=True))

config.sample_len = 16
zoomin_16 = summarize(vp.build_output_filepath(config=config, parent_only=True))


# ZOOM OUT
config.ds_type = 'zoom_out'
config.sample_len = 4
zoomout_4 = summarize(vp.build_output_filepath(config=config, parent_only=True))

config.sample_len = 8
zoomout_8 = summarize(vp.build_output_filepath(config=config, parent_only=True))

config.sample_len = 16
zoomout_16 = summarize(vp.build_output_filepath(config=config, parent_only=True))


# GTSRB
config.ds_type = 'gtsrb'
config.sample_len = 4
gtsrb_4 = summarize(vgp.build_output_filepath(config=config, parent_only=True))

config.sample_len = 8
gtsrb_8 = summarize(vgp.build_output_filepath(config=config, parent_only=True))

config.sample_len = 16
gtsrb_16 = summarize(vgp.build_output_filepath(config=config, parent_only=True))


if __name__ == "__main__":
    plt.figure(figsize=(7, 8))

    plt.plot(x, zoomin_4, linestyle='solid', color='red')
    plt.plot(x, zoomin_8, linestyle='solid', color='green')
    plt.plot(x, zoomin_16, linestyle='solid', color='blue')

    plt.plot(x, zoomout_4, linestyle='dashed', color='red')
    plt.plot(x, zoomout_8, linestyle='dashed', color='green')
    plt.plot(x, zoomout_16, linestyle='dashed', color='blue')

    plt.plot(x, gtsrb_4, linestyle='dotted', color='red')
    plt.plot(x, gtsrb_8, linestyle='dotted', color='green')
    plt.plot(x, gtsrb_16, linestyle='dotted', color='blue')

    plt.yscale('log')

    # legend lines
    line1 = mlines.Line2D([], [], color='black', linestyle='solid', label='Zoom In')
    line2 = mlines.Line2D([], [], color='black', linestyle='dashed', label='Zoom Out')
    line3 = mlines.Line2D([], [], color='black', linestyle='dotted', label='GTSRB')

    # legend colors
    color_red = mpatches.Patch(color='red', label='4 frames')
    color_green = mpatches.Patch(color='green', label='8 frames')
    color_blue = mpatches.Patch(color='blue', label='16 frames')

    # add labels and title
    plt.xlabel("Epsilon")
    plt.ylabel("Runtime (s)")
    plt.title("Comparison of average runtime")
    plt.legend(handles=[line1, line2, line3, color_red, color_green, color_blue], loc='upper left')

    plot_fp = os.path.join(root, 'FORMALISE2025', 'figs', 'runtime_comp.png')
    plt.savefig(plot_fp, dpi=300)
