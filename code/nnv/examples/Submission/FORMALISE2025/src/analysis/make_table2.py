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

# get the results dir
root = os.path.dirname(os.getcwd()) # Submission folder
output_dir = os.path.join(root, 'FORMALISE2025' ,'results')

def summarize(output_file_dir):
    # print(f'{output_file_dir}')
    num_samples_verified = []
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
            num_samples_verified.append(total_verified)

            # calculate average time to verify
            average_time = np.mean(t)

            average_times.append(average_time)

    return num_samples_verified, average_times

if __name__=="__main__":
    # OUTPUTS
    tv = []
    rt = []

    # ZOOM IN RELAX
    config = Config(
            class_size=10,
            epsilon=[1/255, 2/255, 3/255],
            ds_type='zoom_in',
            sample_len=4,
            ver_algorithm='relax',
            timeout=1800,
            output_dir=output_dir
    )
    zoomin_4_tv_relax, zoomin_4_relax = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomin_4_tv_relax)
    rt.append(zoomin_4_relax)

    config.sample_len = 8
    zoomin_8_tv_relax, zoomin_8_relax = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomin_8_tv_relax)
    rt.append(zoomin_8_relax)

    config.sample_len = 16
    zoomin_16_tv_relax, zoomin_16_relax = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomin_16_tv_relax)
    rt.append(zoomin_16_relax)

    # ZOOM IN APPROX
    config.sample_len = 4
    config.ver_algorithm = 'approx'
    zoomin_4_tv_approx, zoomin_4_approx = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomin_4_tv_approx)
    rt.append(zoomin_4_approx)

    config.sample_len = 8
    zoomin_8_tv_approx, zoomin_8_approx = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomin_8_tv_approx)
    rt.append(zoomin_8_approx)


    # ZOOM OUT RELAX
    config.ds_type = 'zoom_out'
    config.sample_len = 4
    config.ver_algorithm = 'relax'
    zoomout_4_tv_relax, zoomout_4_relax = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomout_4_tv_relax)
    rt.append(zoomout_4_relax)

    config.sample_len = 8
    zoomout_8_tv_relax, zoomout_8_relax = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomout_8_tv_relax)
    rt.append(zoomout_8_relax)

    config.sample_len = 16
    zoomout_16_tv_relax, zoomout_16_relax = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomout_16_tv_relax)
    rt.append(zoomout_16_relax)

    # ZOOM OUT APPROX
    config.sample_len = 4
    config.ver_algorithm = 'approx'
    zoomout_4_tv_approx, zoomout_4_approx = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomout_4_tv_approx)
    rt.append(zoomout_4_approx)

    config.sample_len = 8
    zoomout_8_tv_approx, zoomout_8_approx = summarize(vp.build_output_filepath(config=config, parent_only=True))
    tv.append(zoomout_8_tv_approx)
    rt.append(zoomout_8_approx)


    # GTSRB RELAX
    config.ds_type = 'gtsrb'
    config.sample_len = 4
    config.ver_algorithm = 'relax'
    gtsrb_4_tv_relax, gtsrb_4_relax = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(gtsrb_4_tv_relax)
    rt.append(gtsrb_4_relax)

    config.sample_len = 8
    gtsrb_8_tv_relax, gtsrb_8_relax = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(gtsrb_8_tv_relax)
    rt.append(gtsrb_8_relax)

    config.sample_len = 16
    gtsrb_16_tv_relax, gtsrb_16_relax = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(gtsrb_16_tv_relax)
    rt.append(gtsrb_16_relax)

    # GTSRB APPROX
    config.sample_len = 4
    config.ver_algorithm = 'approx'
    gtsrb_4_tv_approx, gtsrb_4_approx = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(gtsrb_4_tv_approx)
    rt.append(gtsrb_4_approx)

    config.sample_len = 8
    gtsrb_8_tv_approx, gtsrb_8_approx = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(gtsrb_8_tv_approx)
    rt.append(gtsrb_8_approx)

    config.sample_len = 16
    gtsrb_16_tv_approx, gtsrb_16_approx = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(gtsrb_16_tv_approx)
    rt.append(gtsrb_16_approx)


    # STMNIST RELAX
    config.ds_type = 'stmnist'
    config.sample_len = 16
    config.ver_algorithm = 'relax'
    stmnist_16_tv_relax, stmnist_16_relax = summarize(vsp.build_output_filepath(config=config, parent_only=True))
    tv.append(stmnist_16_tv_relax)
    rt.append(stmnist_16_relax)

    config.sample_len = 32
    stmnist_32_tv_relax, stmnist_32_relax = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(stmnist_32_tv_relax)
    rt.append(stmnist_32_relax)

    config.sample_len = 64
    stmnist_64_tv_relax, stmnist_64_relax = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(stmnist_64_tv_relax)
    rt.append(stmnist_64_relax)

    # STMNIST APPROX
    config.sample_len = 16
    config.ver_algorithm = 'approx'
    stmnist_16_tv_approx, stmnist_16_approx = summarize(vsp.build_output_filepath(config=config, parent_only=True))
    tv.append(stmnist_16_tv_approx)
    rt.append(stmnist_16_approx)

    config.sample_len = 32
    stmnist_32_tv_approx, stmnist_32_approx = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(stmnist_32_tv_approx)
    rt.append(stmnist_32_approx)

    config.sample_len = 64
    stmnist_64_tv_approx, stmnist_64_approx = summarize(vgp.build_output_filepath(config=config, parent_only=True))
    tv.append(stmnist_64_tv_approx)
    rt.append(stmnist_64_approx)


    if __name__ == "__main__":
        datasets = ["Zoom In", "Zoom Out", "GTSRB", "STMNIST"]    
        veralgs = ["relax", "approx"]

        # Create all the labels
        labels = []
        for ds in datasets:
                for va in veralgs:
                    if (ds == "Zoom In" or ds == "Zoom Out") and va == "approx":
                        sample_lengths = ["4", "8"]

                    elif ds == "STMNIST":
                        sample_lengths = ["4", "8", "16"]
                    else:
                        sample_lengths = ["16", "32", "64"]
                    
                    for sl in sample_lengths:
                        for eps in ["1/255", "2/255", "3/255"]:
                            labels.append(f'{ds} {va} {sl} {eps}')
        
        # needs to be flattened for pairing with label later
        total_verified_flattened = [
            x for sample in tv
            for x in sample
        ]

        # needs to be flattened for pairing with label later
        running_time_flattened = [
            x for sample in rt
            for x in sample
        ]

        MAX_LENGTH_LABEL = max([len(label) for label in labels])

        # Pair together the labels and results
        labels_and_results = zip(labels, total_verified_flattened, running_time_flattened)
        
        # Create the row for the table
        def format_row(sample):
            label, tv, t = sample
            label = f"{label:<{MAX_LENGTH_LABEL}}"
            
            # compute percent verified
            num_attempted = 215 if 'GTSRB' in label else 100
            percent_verified = round((tv / num_attempted) * 100, 2)
            t = round(t, 2)

            pv_string = f"{percent_verified:<{6}}"
            at_string = f"{t:<{9}}"

            return " | ".join([label, pv_string, at_string])
            
        # Create the human-readable table
        output_file = os.path.join(output_dir, "table2.txt")
        with open(output_file, 'w') as f:

            # write the header
            label_header = f"{'LABEL':^{MAX_LENGTH_LABEL}}"
            psrv_header = f"{'PSRV':^{6}}"
            time_header = f"{'Avg. Time':^{10}}"
            header = " | ".join([label_header, psrv_header, time_header])
            line = '-'*len(header)
            f.write(header + "\n")
            f.write(line + "\n")

            for sample in labels_and_results:
                f.write(format_row(sample) + "\n")
                

