import json
import os


def get_averages(filepath):

    epsilon_1_avg = {"correct": 0, "time": 0}
    epsilon_2_avg = {"correct": 0, "time": 0}
    epsilon_3_avg = {"correct": 0, "time": 0}

    with open(filepath, "r") as json_file:
        data = json.load(json_file)

        # iterate over each sample (should be 100)
        for i in range(100):
            res = data["res"][i]
            t = data["time"][i]

            # Certified robustness accuracy
            epsilon_1_avg["correct"] += res[0] if res[0] == 1 else 0
            epsilon_2_avg["correct"] += res[1] if res[1] == 1 else 0
            epsilon_3_avg["correct"] += res[2] if res[2] == 1 else 0

            # Compute time
            epsilon_1_avg["time"] += t[0]
            epsilon_2_avg["time"] += t[1]
            epsilon_3_avg["time"] += t[2]

        # compute CRA
        cra = {
            "epsilon_1": epsilon_1_avg["correct"] / 100,
            "epsilon_2": epsilon_2_avg["correct"] / 100,
            "epsilon_3": epsilon_3_avg["correct"] / 100,
        }

        # compute avg comp time
        compute_time = {
            "epsilon_1": epsilon_1_avg["time"] / 100,
            "epsilon_2": epsilon_2_avg["time"] / 100,
            "epsilon_3": epsilon_3_avg["time"] / 100,
        }

    return cra, compute_time


def get_scalability_averages(filepath):

    epsilon_1_avg = {"correct": 0, "time": 0}
    epsilon_2_avg = {"correct": 0, "time": 0}
    epsilon_3_avg = {"correct": 0, "time": 0}

    epsilon_1_per_frame = [{"correct": 0, "time": 0} for _ in range(8)]
    epsilon_2_per_frame = [{"correct": 0, "time": 0} for _ in range(8)]
    epsilon_3_per_frame = [{"correct": 0, "time": 0} for _ in range(8)]

    with open(filepath, "r") as json_file:
        data = json.load(json_file)

        # iterate over each sample (should be 10)
        for i in range(10):
            res = data["res"][i]
            t = data["time"][i]

            # Certified robustness accuracy
            epsilon_1_avg["correct"] += 1 if all(res[0]) else 0
            epsilon_2_avg["correct"] += 1 if all(res[1]) else 0
            epsilon_3_avg["correct"] += 1 if all(res[2]) else 0

            # Compute time
            epsilon_1_avg["time"] += sum(t[0])
            epsilon_2_avg["time"] += sum(t[1])
            epsilon_3_avg["time"] += sum(t[2])

            # Per frame metrics
            for frame_i in range(len(res[0])):  # epsilon 1
                epsilon_1_per_frame[frame_i]["correct"] += (
                    1 if res[0][frame_i] == 1 else 0
                )
                epsilon_1_per_frame[frame_i]["time"] += t[0][frame_i]

            for frame_i in range(len(res[1])):  # epsilon 2
                epsilon_2_per_frame[frame_i]["correct"] += (
                    1 if res[1][frame_i] == 1 else 0
                )
                epsilon_2_per_frame[frame_i]["time"] += t[1][frame_i]

            for frame_i in range(len(res[2])):  # epsilon 3
                epsilon_3_per_frame[frame_i]["correct"] += (
                    1 if res[2][frame_i] == 1 else 0
                )
                epsilon_3_per_frame[frame_i]["time"] += t[2][frame_i]

        # compute CRA
        cra = {
            "epsilon_1": epsilon_1_avg["correct"] / 10,
            "epsilon_2": epsilon_2_avg["correct"] / 10,
            "epsilon_3": epsilon_3_avg["correct"] / 10,
        }

        # compute avg comp time
        compute_time = {
            "epsilon_1": epsilon_1_avg["time"] / 10,
            "epsilon_2": epsilon_2_avg["time"] / 10,
            "epsilon_3": epsilon_3_avg["time"] / 10,
        }

        # compute per frame CRA
        per_frame_cra = {
            "epsilon_1": [
                epsilon_1_per_frame[frame_i]["correct"] / 10
                for frame_i in range(len(res[0]))
            ],
            "epsilon_2": [
                epsilon_2_per_frame[frame_i]["correct"] / 10
                for frame_i in range(len(res[1]))
            ],
            "epsilon_3": [
                epsilon_3_per_frame[frame_i]["correct"] / 10
                for frame_i in range(len(res[2]))
            ],
        }

        # compute per frame avg comp time
        per_frame_time = {
            "epsilon_1": [
                epsilon_1_per_frame[frame_i]["time"] / 10
                for frame_i in range(len(res[0]))
            ],
            "epsilon_2": [
                epsilon_2_per_frame[frame_i]["time"] / 10
                for frame_i in range(len(res[1]))
            ],
            "epsilon_3": [
                epsilon_3_per_frame[frame_i]["time"] / 10
                for frame_i in range(len(res[2]))
            ],
        }

    return cra, compute_time, per_frame_cra, per_frame_time


if __name__ == "__main__":

    filepaths = [
        "src/vvn/results/ZoomOut/C3D_small_robustness_results_FPV.json",
        "src/vvn/results/ZoomIn/C3D_small_robustness_results_FPV.json",
        "src/vvn/results/ScalabilityZoomOut/C3D_small_robustness_results_FPV.json",
        "src/vvn/results/ScalabilityZoomIn/C3D_small_robustness_results_FPV.json",
    ]

    results = {}

    print(os.getcwd())

    for fp in filepaths:
        shortName = fp.split("/")[3]
        print("Getting averages for ... " + shortName)

        # first two are normal averages
        if "Scalability" not in fp:
            cra, t = get_averages(fp)
            results[shortName] = {"cra": cra, "time": t}
        else:
            cra, t, per_frame_cra, per_frame_t = get_scalability_averages(fp)
            results[shortName] = {
                "cra": cra,
                "time": t,
                "per_frame_cra": per_frame_cra,
                "per_frame_time": per_frame_t,
            }

    with open("src/analysis/final_results.json", "w") as final:
        json.dump(results, final, indent=4)

    print("All done!")
