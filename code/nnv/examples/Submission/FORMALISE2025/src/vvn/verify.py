# python stdlib 
import csv
import io
import os
import sys
from typing import Tuple

# third-party packages
from dotenv import load_dotenv
import matlab.engine
import numpy as np

# local modules
import vvn.prep as vp
from vvn.config import Config
import vvn.gtsrbprep as vgp
import vvn.stmnistprep as vsp

# for SoSym (on docker)
RESULTS_DIR = '/tmp'

# load environment variables 
load_dotenv()

# define global variables
NNV_PATH = os.environ['NNV_PATH']
NPY_MATLAB_PATH = os.environ['NPY_MATLAB_PATH']

def prepare_engine(nnv_path, npy_matlab_path):
    if not nnv_path or not npy_matlab_path:
        raise Exception('One of nnv_path or npy_matlab_path is not defined. Please ensure these have been set before running.')

    # start matlab
    eng = matlab.engine.start_matlab()
    print('started matlab engine!')

    # add nnv path and npy-matlab path
    eng.addpath(os.getcwd())
    eng.addpath(eng.genpath(nnv_path))
    eng.addpath(eng.genpath(npy_matlab_path))

    # save reference to it for calling matlab scripts to engine later
    return eng

def verify(ds_type, sample_len, ver_algorithm, eng, index, eps_index, timeout) -> Tuple[int, float | str, str]:
    # check that MATLAB engine was started correctly and is accessible
    if not eng:
        raise Exception('MATLAB Engine was not correctly started and shared. Please make sure to run `prepare_engine`.')

    # call to MATLAB script to run verification
    future = eng.verifyvideo(ds_type, sample_len, ver_algorithm, index, eps_index, nargout=3, background=True, stdout=io.StringIO())

    try:
        [res, t, met] = future.result(timeout=float(timeout))

    except matlab.engine.TimeoutError:
        print('timeout')
        res = 3
        t = 'timeout'
        met = 'timeout'

    future.cancel()

    return res, t, met

def run(config, indices) -> None:
    # Unpack configuration settings;
    epsilon = config.epsilon
    timeout = config.timeout
    
    ds_type = config.ds_type
    sample_len = config.sample_len
    ver_algorithm = config.ver_algorithm
    output_dir = config.output_dir

    print(f'Running verification with config: verification algorithm={ver_algorithm}, dataset type={ds_type}, video length={sample_len}') 

    # make sure directory structure + results files are created and correct
    # vp.prepare_filetree(config)

    # make sure matlab is started
    eng = prepare_engine(NNV_PATH, NPY_MATLAB_PATH)

    # start verification
    for sample_num, index in enumerate(indices):
        print(f'Iteration {sample_num + 1}')

        # select epsilon
        for eps_index in range(1, len(epsilon) + 1):
            # output_file = vp.build_output_filepath(config=config, filename=f'eps={eps_index}_255')
            output_file_ds_type = 'ZoomIn' if ds_type == 'zoom_in' else 'ZoomOut'
            output_file = os.path.join(output_dir, output_file_ds_type, str(sample_len), f'eps={eps_index}_255.csv')

            # verify the sample with a specific epsilon value
            res, t, met = verify(ds_type, sample_len, ver_algorithm, eng, index, eps_index, timeout)

            # write the results
            write_results(output_file, sample_num, res, t, met)

    # close matlab after experiment finishes
    eng.quit()

def verify_gtsrb(sample_len, ver_algorithm, eng, index, eps_index, timeout) -> Tuple[int, float | str, str]:
    # check that MATLAB engine was started correctly and is accessible
    if not eng:
        raise Exception('MATLAB Engine was not correctly started and shared. Please make sure to run `prepare_engine`.')

    # call to MATLAB script to run verification
    future = eng.verifygtsrb(sample_len, ver_algorithm, index, eps_index, nargout=3, background=True, stdout=io.StringIO())

    try:
        [res, t, met] = future.result(timeout=float(timeout))

    except matlab.engine.TimeoutError:
        print('timeout')
        res = 3
        t = 'timeout'
        met = 'timeout'

    future.cancel()

    return res, t, met

def run_gtsrb(config, indices) -> None:
    # Unpack configuration settings;
    epsilon = config.epsilon
    timeout = config.timeout
    
    ds_type = config.ds_type
    sample_len = config.sample_len
    ver_algorithm = config.ver_algorithm
    output_dir = config.output_dir

    print(f'Running verification with config: verification algorithm={ver_algorithm}, dataset type={ds_type}, video length={sample_len}') 

    # make sure directory structure + results files are created and correct
    # vgp.prepare_filetree(config)

    # make sure matlab is started
    eng = prepare_engine(NNV_PATH, NPY_MATLAB_PATH)

    # start verification
    for sample_num, index in enumerate(indices):
        print(f'Iteration {sample_num + 1}')

        if_timeout = False

        # select epsilon
        for eps_index in range(1, len(epsilon) + 1):
            # output_file = vgp.build_output_filepath(config=config, filename=f'eps={eps_index}_255')
            output_file = os.path.join(output_dir, str(sample_len), f'eps={eps_index}_255.csv')

            # skip if timeout was met at any point in the previous iterations
            if if_timeout:
                res, t, met = 3, "timeout", "timeout"
                write_results(output_file, sample_num, res, t, met)
                continue

            # verify the sample with a specific epsilon value
            res, t, met = verify_gtsrb(sample_len, ver_algorithm, eng, index, eps_index, timeout)

            if res == 3:
                if_timeout = True

            # write the results
            write_results(output_file, sample_num, res, t, met)

    # close matlab after experiment finishes
    eng.quit()

def verify_stmnist(sample_len, ver_algorithm, eng, index, eps_index, timeout) -> Tuple[int, float | str, str]:
    # check that MATLAB engine was started correctly and is accessible
    if not eng:
        raise Exception('MATLAB Engine was not correctly started and shared. Please make sure to run `prepare_engine`.')

    # call to MATLAB script to run verification
    future = eng.verifystmnist(sample_len, ver_algorithm, index, eps_index, nargout=3, background=True, stdout=io.StringIO())

    try:
        [res, t, met] = future.result(timeout=float(timeout))

    except matlab.engine.TimeoutError:
        print('timeout')
        res = 3
        t = 'timeout'
        met = 'timeout'

    future.cancel()

    return res, t, met

def run_stmnist(config, indices) -> None:
    # Unpack configuration settings;
    epsilon = config.epsilon
    timeout = config.timeout
    
    ds_type = config.ds_type
    sample_len = config.sample_len
    ver_algorithm = config.ver_algorithm
    output_dir = config.output_dir

    print(f'Running verification with config: verification algorithm={ver_algorithm}, dataset type={ds_type}, video length={sample_len}') 

    # make sure directory structure + results files are created and correct
    # vsp.prepare_filetree(config)

    # make sure matlab is started
    eng = prepare_engine(NNV_PATH, NPY_MATLAB_PATH)

    # start verification
    for sample_num, index in enumerate(indices):
        print(f'Iteration {sample_num + 1}')

        if_timeout = False

        # select epsilon
        for eps_index in range(1, len(epsilon) + 1):
            # output_file = vsp.build_output_filepath(config=config, filename=f'eps={eps_index}_255')
            output_file = os.path.join(output_dir, str(sample_len), f'eps={eps_index}_255.csv')

            # skip if timeout was met at any point in the previous iterations
            if if_timeout:
                res, t, met = 3, "timeout", "timeout"
                write_results(output_file, sample_num, res, t, met)
                continue

            # verify the sample with a specific epsilon value
            res, t, met = verify_stmnist(sample_len, ver_algorithm, eng, index, eps_index, timeout)

            if res == 3:
                if_timeout = True

            # write the results
            write_results(output_file, sample_num, res, t, met)

    # close matlab after experiment finishes
    eng.quit()

def verify_kthactions(sample_len, ver_algorithm, eng, index, eps_index, timeout) -> Tuple[int, float | str, str]:
    # check that MATLAB engine was started correctly and is accessible
    if not eng:
        raise Exception('MATLAB Engine was not correctly started and shared. Please make sure to run `prepare_engine`.')

    # call to MATLAB script to run verification
    future = eng.verifykthactions(sample_len, ver_algorithm, index, eps_index, nargout=3, background=True, stdout=io.StringIO())

    try:
        [res, t, met] = future.result(timeout=float(timeout))

    except matlab.engine.TimeoutError:
        print('timeout')
        res = 3
        t = 'timeout'
        met = 'timeout'

    future.cancel()

    return res, t, met

def run_kthactions(config, indices) -> None:
    # Unpack configuration settings;
    epsilon = config.epsilon
    timeout = config.timeout
    
    ds_type = config.ds_type
    sample_len = config.sample_len
    ver_algorithm = config.ver_algorithm

    print(f'Running verification with config: verification algorithm={ver_algorithm}, dataset type={ds_type}, video length={sample_len}') 

    # make sure directory structure + results files are created and correct
    results_dir = os.path.join(RESULTS_DIR, 'KTHActions')

    if not os.path.isdir(results_dir):
        os.mkdir(results_dir)

    # make results dir for videos of specific sample length
    if not os.path.isdir(os.path.join(results_dir, str(sample_len))):
        os.mkdir(os.path.join(results_dir, str(sample_len)))

    # make sure matlab is started
    eng = prepare_engine(NNV_PATH, NPY_MATLAB_PATH)

    # start verification
    for sample_num, index in enumerate(indices):
        print(f'Iteration {sample_num + 1}')

        if_timeout = False

        # select epsilon
        for eps_index in range(1, len(epsilon) + 1):
            output_file = os.path.join(results_dir, str(sample_len), f'eps={eps_index}_255.csv')

            # skip if timeout was met at any point in the previous iterations
            if if_timeout:
                res, t, met = 3, "timeout", "timeout"
                write_results(output_file, sample_num, res, t, met)
                continue

            # verify the sample with a specific epsilon value
            res, t, met = verify_kthactions(sample_len, ver_algorithm, eng, index, eps_index, timeout)

            if res == 3:
                if_timeout = True

            # write the results
            write_results(output_file, sample_num, res, t, met)

    # close matlab after experiment finishes
    eng.quit()

def verify_ucf11(sample_len, ver_algorithm, eng, index, eps_index, timeout, out_channels) -> Tuple[int, float | str, str]:
    # check that MATLAB engine was started correctly and is accessible
    if not eng:
        raise Exception('MATLAB Engine was not correctly started and shared. Please make sure to run `prepare_engine`.')

    # call to MATLAB script to run verification
    # future = eng.verifyucf11(sample_len, ver_algorithm, index, eps_index, nargout=3, background=True, stdout=io.StringIO())
    buf = io.StringIO()
    future = eng.verifyucf11(sample_len, ver_algorithm, index, eps_index, out_channels, nargout=3, background=True, stdout=buf, stderr=buf)
    print(buf.getvalue())

    try:
        [res, t, met] = future.result(timeout=float(timeout))

    except matlab.engine.TimeoutError:
        print('timeout')
        res = 3
        t = 'timeout'
        met = 'timeout'

    future.cancel()

    return res, t, met

def run_ucf11(config, indices, out_channels) -> None:
    # Unpack configuration settings;
    epsilon = config.epsilon
    timeout = config.timeout
    
    ds_type = config.ds_type
    sample_len = config.sample_len
    ver_algorithm = config.ver_algorithm

    print(f'Running verification with config: verification algorithm={ver_algorithm}, dataset type={ds_type}, video length={sample_len}') 

    # make sure directory structure + results files are created and correct
    results_dir = os.path.join(RESULTS_DIR, 'UCF11')

    if not os.path.isdir(results_dir):
        os.mkdir(results_dir)

    # make results dir for videos of specific sample length if it doesn't exist
    if not os.path.isdir(os.path.join(results_dir, str(sample_len))):
        os.mkdir(os.path.join(results_dir, str(sample_len)))

    # make results dir for videos with models of specific # of out channels if it doesn't exist
    if not os.path.isdir(os.path.join(results_dir, str(sample_len), f'{out_channels}outchannels')):
        os.mkdir(os.path.join(results_dir, str(sample_len), f'{out_channels}outchannels'))

    # make sure matlab is started
    eng = prepare_engine(NNV_PATH, NPY_MATLAB_PATH)

    # start verification
    for sample_num, index in enumerate(indices):
        print(f'Iteration {sample_num + 1}')

        if_timeout = False

        # select epsilon
        for eps_index in range(1, len(epsilon) + 1):
            output_file = os.path.join(results_dir, str(sample_len), f'{out_channels}outchannels', f'eps={eps_index}_255.csv')

            # skip if timeout was met at any point in the previous iterations
            if if_timeout:
                res, t, met = 3, "timeout", "timeout"
                write_results(output_file, sample_num, res, t, met)
                continue

            # verify the sample with a specific epsilon value
            res, t, met = verify_ucf11(sample_len, ver_algorithm, eng, index, eps_index, timeout, out_channels)

            if res == 3:
                if_timeout = True

            # write the results
            write_results(output_file, sample_num, res, t, met)

    # close matlab after experiment finishes
    eng.quit()

def write_results(output_file, sample_num, res, t, met):
    with open(output_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([sample_num, res, t, met])

def summarize(output_file_dir, data_len):
    print(f'{output_file_dir}')
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

        # display the results
        results_header_str = f'Results of verification with {filename.split(".")[0]}'
        total_verified_str = f'Verified {int(total_verified)} robust samples out of {data_len}.'
        average_time_str = f'Average running time was : {average_time}.'
        rowlength = max(len(total_verified_str), len(average_time_str), len(results_header_str))
        print('='*rowlength)
        print(results_header_str)
        print('---')
        print(total_verified_str)
        print('---')
        print(average_time_str)
        print('='*rowlength)
        print('\n\n')


if __name__ == "__main__":
    # for smoke test
    dst = sys.argv[1]
    veralg = sys.argv[2]
    epsilon_index = sys.argv[3]
    sample_len = sys.argv[4]
    timeout = sys.argv[5]
    index = sys.argv[6]

    # verify command line arguments
    if dst not in ['zoom_in', 'zoom_out', 'gtsrb', 'stmnist']:
        raise error("Invalid dataset argument. Must be one of 'zoom_in', 'zoom_out', 'gtsrb' or 'stmnist'.")
    
    if veralg not in ['relax', 'approx']:
        raise error("Invalid verification algorithm. Must be one of 'relax' or 'approx'.")
    
    if epsilon_index not in ['1', '2', '3']:
        raise error("Invalid epsilon. Must be one of 1, 2, or 3.")

    if dst == 'stmnist':
        if sample_len not in ['16', '32', '64']:
            raise error("Invalid sample length for STMNIST. Must be one of 16, 32, 64.")
    else:
        if sample_len not in ['4', '8', '16']:
            raise error("Invalid sample length for MNIST and GTSRB Video. Must be one of 4, 8, 16.")

    # convert CL arguments to correct data types
    epsilon_index = int(epsilon_index)
    sample_len = int(sample_len)
    index = int(index)
    timeout = int(timeout)

    # start matlab
    eng = prepare_engine(NNV_PATH, NPY_MATLAB_PATH)

    if dst == 'zoom_in' or dst == 'zoom_out':
        res, t, met = verify(dst, sample_len, veralg, eng, index, epsilon_index, timeout)
    elif dst == 'stmnist':
        res, t, met = verify_stmnist(sample_len, veralg, eng, index, epsilon_index, timeout)
    else:
        res, t, met = verify_gtsrb(sample_len, veralg, eng, index, epsilon_index, timeout)

    print(f'{dst}-{sample_len}f with {veralg} and e={epsilon_index}/255.') 
    print(f'Res: {res}, Time: {t}, Met: {met}\n')
