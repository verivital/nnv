# sys.path.insert(0, BASE_FOLDER)
import os
import pickle
import sys
from collections import defaultdict
from datetime import datetime
from functools import partial
from timeit import default_timer as timer
from typing import Tuple, List, Union

import numpy as np
from prettytable import PrettyTable
from tqdm import tqdm

from RNN.Adversarial import adversarial_query, get_out_idx
from RNN.MultiLayerBase import GurobiMultiLayer
from rnn_experiment.self_compare.create_sbatch_iterations_exp import BASE_FOLDER, OUT_FOLDER

MODELS_FOLDER = os.path.join(BASE_FOLDER, "models/")
POINTS_PATH = os.path.join(MODELS_FOLDER, "points.pkl") 

IN_SHAPE = (40,)
NUM_SAMPLE_POINTS = 25
NUM_RUNNER_UP = 1
ALL_NETS = [f for f in os.listdir(MODELS_FOLDER) if f.endswith(".h5")]


def run_experiment(in_tensor, radius, idx_max, other_idx, h5_file, gurobi_ptr, n_iterations, steps):
    queries_stats = {}
    start = timer()
    try:
        res, queries_stats, alpha_history = adversarial_query(in_tensor, radius, idx_max, other_idx, h5_file,
                                                              gurobi_ptr, n_iterations, steps)
    except ValueError as e:
        # row_result = {'point': in_tensor, 'error': e, 'error_traceback': traceback.format_exc(), 'result' : False}
        res = False
    except TimeoutError as e:
        res = False
        queries_stats['FFNN_Timeout'] = True
    except AssertionError as e:
        res = False

    end = timer()

    if queries_stats is not None:
        queries_stats['total_time'] = end - start
    return {'time': end - start, 'result': res, 'stats': queries_stats}


def run_all_experiments(net_options, points, t_range, other_idx_method, gurobi_ptr, radius=0.10, steps_num=1500,
                        save_results=True):
    print('================ eps: %s ================', radius)
    # assert len(points) > 20
    results = defaultdict(list)
    if len(net_options) == 1:
        net_name = ''.join(net_options[0].split('.')[:-1]).split('/')[-1]
    else:
        net_name = ''
    pickle_path = os.path.join(OUT_FOLDER,
                               'gurobi' + str(datetime.now()).replace('.', '').replace(' ', '') + "{}.pkl".format(
                                   net_name))
    print("starting fresh experiment", "\n", "#" * 100)
    partial_results = {}

    print("#" * 100, "\nwriting results to: {}".format(pickle_path), "\n", "#" * 100)
    counter = 0
    pbar = tqdm(total=len(other_idx_method) * len(points) * len(net_options) * len(t_range))
    for method in other_idx_method:
        for idx, point in enumerate(points):
            for path in net_options:
                if not os.path.exists(path):
                    path = os.path.join(MODELS_FOLDER, path)
                    if not os.path.exists(path):
                        raise FileNotFoundError(path)
                for t in t_range:
                    if counter < 0:
                        counter += 1
                        pbar.update(1)
                        have_point = True
                    else:
                        have_point = False
                        net_name = ''.join(path.split('.')[:-1]).split('/')[-1]
                        name = "{}_{}_{}".format(net_name, radius, t)

                        if name in partial_results:
                            for res in partial_results[name]:
                                if not have_point and res['t'] == t and \
                                        (('in_tensor' in res and np.all(res['in_tensor'] == point)) or
                                         ('in_tesnor' in res and np.all(res['in_tesnor'] == point))):
                                    # already have this result
                                    pbar.update(1)
                                    results[name].append(res)
                                    have_point = True
                    if not have_point:
                        idx_max, other_idx = get_out_idx(point, t, path, method)
                        net_name = ''.join(path.split('.')[:-1]).split('/')[-1]

                        result = run_experiment(point, radius, idx_max, other_idx, path, gurobi_ptr, t, steps=steps_num)
                        result.update({'h5_file': net_name, 't': t, 'other_idx': other_idx, 'in_tensor': point,
                                       'steps_num': steps_num})
                        results[name].append(result)
                        if not result['result']:
                            print("FAIL on point index: {}".format(idx))
                        pbar.update(1)
                        if save_results:
                            pickle.dump(results, open(pickle_path, "wb"))
    if save_results:
        parse_results_file(pickle_path, radius, print_latex=True)
    return results


def extract_model_name(name: str) -> str:
    if 'rnn' not in name or 'fc32' not in name:
        return name
    model_name = 'rnn'
    fc_count = 0
    for w in name.split('_'):
        if 'rnn' in w:
            model_name += w.replace('rnn', '')
        elif 'fc' in w:
            fc_count += 1
        else:
            model_name += '_{}fc32'.format(fc_count)
    return model_name


def parse_results_file(name_path_map: Union[List[Tuple[str, str]], str], radius, print_latex=False):
    if isinstance(name_path_map, str) and name_path_map.endswith('pkl'):
        os.makedirs("temp/", exist_ok=True)
        results = pickle.load(open(name_path_map, "rb"))
        name_path_map = []
        new_files = defaultdict(dict)
        for k, v in results.items():
            model_name = '_'.join(k.split('_')[:-1])
            new_files[model_name].update({k: v})
        for k, v in new_files.items():
            pickle.dump(v, open("temp/{}".format(k), "wb"))
            name_path_map.append((k.replace("model_20classes_", ""), "temp/{}".format(k)))
    if isinstance(name_path_map, str) and os.path.isdir(name_path_map):
        tuples = []
        for f in sorted(os.listdir(name_path_map)):
            if not f.endswith('.pkl'):
                continue
            p = os.path.join(name_path_map, f)
            if not os.path.isfile(p):
                continue
            m = extract_model_name(f)
            tuples.append((m, p))
        name_path_map = tuples

    x = PrettyTable()
    if len(set(name_path_map)) != len(name_path_map):
        print('PROBLEM')
    x.field_names = ['Tmax'] + [extract_model_name(p[0]) for p in name_path_map]  # if p[0] not in names]
    rows = {}

    total_time = 0
    total_points = 0
    total_timeout = 0
    total_gurobi_time = 0
    total_invariant_time = 0
    total_property_time = 0
    total_init_time = 0
    sum_success = 0
    max_query_time = -1
    all_times = []
    for p in name_path_map:
        d = pickle.load(open(p[1], "rb"))
        for key, value in d.items():
            net_name = "_".join(key.split("_")[:-2])
            t = int(key.split("_")[-1])
            assert t >= 2
            if t not in rows:
                rows[t] = [t]

            res = parse_dictionary(value)
            success_rate = int(res['success_rate'] * 100)
            all_times += res['time']
            if max(res['time']) > max_query_time:
                max_query_time = max(res['time'])
            total_success = res['total_success']
            sum_success += total_success
            avg_run_time = res['avg_total_time_no_timeout']
            avg_init_time = res['avg_initialize_query_time']

            total_time += avg_run_time * res['total']
            # assert np.abs(total_time - res['avg_total_time']) < 10 ** -3
            total_points += res['total']
            total_timeout += res['len_timeout']
            total_init_time += res['avg_initialize_query_time'] * res['total']
            total_gurobi_time += res['avg_step_time'] * res['total']
            total_invariant_time += res['avg_invariant_time'] * res['total']
            total_property_time += res['avg_property_time'] * res['total']
            # assert timeout == 0
            gurobi_time = res['avg_step_time_no_timeout']
            ffnn_time = res['avg_invariant_time'] + res['avg_property_time'] + gurobi_time
            assert ffnn_time < avg_run_time, "{}, {}, {}".format(ffnn_time, gurobi_time, avg_run_time)
            # print(p)s
            # print("total time: {}".format(total_time ))
            # print("gurobi time: {}".format(total_gurobi_time))
            # print("prop time: {}".format(total_property_time))
            # print("ffnn time: {}".format(total_invariant_time + total_property_time))
            # avg_marabou_invariant = res['avg_invariant_time_no_timeout'] - avg_gurobi_invariant
            # print("Format is: time (#success/#total) (#timeout)")
            # rows[t - t_range[0]].append("%.2f (%.2f,%.2f) %d/%d (%d)" % (avg_run_time, ffnn_time, gurobi_time,
            #                                                              total_success, res['total'],
            #                                                              res['len_timeout']))
            rows[t].append("%.2f(%d/%d)" % (avg_run_time - avg_init_time, total_success, res['total']))

    for row in rows.values():
        x.add_row(row)
    # for k, v in rows.items():
    #     x.add_row([k] + ["{}({})".format(v[net_name][0], v[net_name][1]) for net_name in x.field_names])
    print("Format is: time (#success/#total)")
    print(x)
    if print_latex:
        for t, row in enumerate(rows.values()):
            print("\t{} &&&".format(t + 2))
            print(" &&& ".join(row[1:]).replace("  ", " "), "&")
            print("\t\\\\")

    avg_time = total_time / total_points
    avg_gurobi = total_gurobi_time / total_points
    avg_marabou_invariant = total_invariant_time / total_points
    avg_property = total_property_time / total_points
    avg_init = total_init_time / total_points
    print("#" * 100)

    avg_actual_time = avg_time - avg_init
    print("Average run time {:.2f} seconds ({:.2f} including creating the query), over {} points, timeout: {}, success: {} ({:.2f})"
          .format(avg_actual_time, avg_time, total_points, total_timeout, sum_success, sum_success / total_points))
    print("Max query took: {:.2f} seconds ({:.2f} minutes), median: {:.2f}, 90percentie: {:.2f}".format(
        max_query_time, max_query_time / 60, np.median(all_times), np.percentile(all_times, 90)))
    print("avg Time in Gurobi: {:.2f}({:.2f}%), avg Time proving Invariant in Marabou {:.2f} ({:.2f}%),"
          "avg Time proving property: {:.2f} ({:.2f}%), avg initialize time: {:.2f}"
          .format(avg_gurobi, avg_gurobi / avg_actual_time,
                  avg_marabou_invariant, avg_marabou_invariant / avg_actual_time,
                  avg_property, avg_property / avg_actual_time,
                  avg_init))
    print("radius / epsilon: ", radius)
    print("#" * 100)


def parse_dictionary(exp):
    # This is an experiment entry strcture:
    #    { 'time', 'result', 'h5_file', 't', other_idx', 'in_tesnor',  'steps_num'
    #   'stats':
    #       {'property_times': {'avg', 'median', 'raw'},
    #       'invariant_times': {'avg', 'median', 'raw'},
    #       'step_times': {'avg', 'median', 'raw'},
    #       'step_querites', 'property_queries', 'invariant_queries', 'number_of_updates', 'total_time'},
    #   }

    success_exp = [e for e in exp if e['result']]
    timeout_exp = []
    no_timeout_exp = []
    for e in exp:
        if 'number_of_updates' in e['stats'] and e['stats']['number_of_updates'] == e['stats']['property_queries'] \
                and not e['result']:
            timeout_exp.append(e)
        elif 'FFNN_Timeout' in e['stats']:
            timeout_exp.append(e)
        else:
            no_timeout_exp.append(e)

    safe_mean = lambda x: np.mean(x) if len(x) > 0 else 0

    d = {
        'total': len(exp),
        'total_success': len(success_exp),
        'success_rate': len(success_exp) / len(exp),
        'len_timeout': len(timeout_exp),
        'time': [e['time'] - e['stats']['query_initialize'] for e in exp],
        'invariant_time': [sum(e['stats']['invariant_times']['raw']) for e in exp],
        'gurobi_time': [sum(e['stats']['step_times']['raw']) for e in exp],
        'avg_total_time': safe_mean([e['time'] for e in exp]),
        'avg_total_time_no_timeout': safe_mean([e['time'] for e in no_timeout_exp]),
        'avg_invariant_time_no_timeout': safe_mean(
            [sum(e['stats']['invariant_times']['raw']) for e in no_timeout_exp if 'invariant_times' in e['stats']]),
        'avg_property_time_no_timeout': safe_mean(
            [e['stats']['property_times']['avg'] for e in no_timeout_exp if 'property_times' in e['stats']]),
        'avg_step_time_no_timeout': safe_mean(
            [e['stats']['step_times']['avg'] for e in no_timeout_exp if 'step_times' in e['stats']]),
        'avg_invariant_time': safe_mean(
            [sum(e['stats']['invariant_times']['raw']) for e in exp if 'invariant_times' in e['stats']]),
        'avg_property_time': safe_mean(
            [sum(e['stats']['property_times']['raw']) for e in exp if 'property_times' in e['stats']]),
        'avg_step_time': safe_mean([sum(e['stats']['step_times']['raw']) for e in exp if 'step_times' in e['stats']]),
        'num_invariant_avg': safe_mean(
            [e['stats']['invariant_queries'] for e in exp if 'invariant_queries' in e['stats']]),
        'num_property_avg': safe_mean(
            [e['stats']['property_queries'] for e in exp if 'property_queries' in e['stats']]),
        'num_step_avg': safe_mean([e['stats']['step_queries'] for e in exp if 'step_queries' in e['stats']]),
        'avg_total_time_success': safe_mean([e['time'] for e in success_exp]),
        'avg_invariant_time_success': safe_mean([sum(e['stats']['invariant_times']['raw']) for e in success_exp]),
        'avg_property_time_success': safe_mean([sum(e['stats']['property_times']['raw']) for e in success_exp]),
        'avg_step_time_success': safe_mean([sum(e['stats']['step_times']['raw']) for e in success_exp]),
        'num_invariant_avg_success': safe_mean([e['stats']['invariant_queries'] for e in success_exp]),
        'num_property_avg_success': safe_mean([e['stats']['property_queries'] for e in success_exp]),
        'num_step_avg_success': safe_mean([e['stats']['step_queries'] for e in success_exp]),
        'avg_initialize_query_time': safe_mean(
            [e['stats']['query_initialize'] for e in success_exp if 'query_initialize' in e['stats']])
    }

    gurobi_time = d['avg_step_time_no_timeout']
    ffnn_time = d['avg_invariant_time_no_timeout'] + d['avg_property_time_no_timeout'] - gurobi_time
    avg_run_time = d['avg_total_time_no_timeout']

    assert ffnn_time < avg_run_time, "{}, {}, {}".format(ffnn_time, gurobi_time, avg_run_time)
    return d


def parse_inputs(points, other_idx_method, gurobi_ptr):
    save_results = True

    if sys.argv[1] == 'analyze':
        if len(sys.argv) < 2:
            print("USAGE analyze <path> --> where path is a path to single pickle results file or directory with multiple")
        parse_results_file(sys.argv[2], print_latex=0)
    elif sys.argv[1] == 'exp':
        if len(sys.argv) > 3:
            if sys.argv[2] == 'all':
                net = ALL_NETS
            else:
                net = [sys.argv[2]]
            t_range = range(2, int(sys.argv[3]) + 1)  # make T_max inclusive
            run_all_experiments(net, points, t_range, other_idx_method, gurobi_ptr, steps_num=10,
                                save_results=save_results)
        else:
            print("USAGE: exp <NET_PATH or all > <T_max>")
            exit(1)
    else:
        print("USAGE: <exp  or analyze > followed by net + T_max or path respectively") 
        exit(1)


if __name__ == "__main__":
    other_idx_method = [lambda x: np.argsort(x)[-i] for i in range(2, 2 + NUM_RUNNER_UP)]

    points = pickle.load(open(POINTS_PATH, "rb"))[:NUM_SAMPLE_POINTS]
    gurobi_ptr = partial(GurobiMultiLayer, polyhedron_max_dim=2, use_relu=True, add_alpha_constraint=True,
                         use_counter_example=True, debug=1)

    if len(sys.argv) > 1:
        net_options = None
        parse_inputs(points, other_idx_method, gurobi_ptr)
        exit(0)

    t_range = [5]
    radius = [0.02, 0.03, 0.04, 0.05]
    nets = [f for f in os.listdir(MODELS_FOLDER) if "20classes_rnn8" in f]

    for r in radius:
        print('================ eps: %s ================1', r)
        run_all_experiments(nets, points, t_range, other_idx_method, gurobi_ptr, radius = r, save_results=True, steps_num=1)
    exit(0)


