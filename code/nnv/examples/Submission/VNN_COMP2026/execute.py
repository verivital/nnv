"""Launch MATLAB Engines with Python

This module was created for local testing of reachability analysis (https://github.com/verivital/nnv.git) ahead of VNN-COMP 2025.
To perform testing, the following functionality is enabled.

- Start a MATLAB engine.
- Run an instance of reachability analysis given some arguments through the open MATLAB engine.

Arguments can be passed to the call to the script and will be parsed accordingly.
"""


import sys
import matlab.engine
import time
import os


def prepare_instance(category: str, onnx: str, vnnlib: str) -> None:
    """Set up the MATLAB engine for running an instance.

    Parameters:
        onnx   (str): the path to the .onnx file
        vnnlib (str): the path to the .vnnlib file
    """

    print("We should not be here...")


def run_instance(category, onnx, vnnlib, timeout, outputlocation) -> None:
    """Run an instance based on parameters defined in .csv file.

    Parameters:
        onnx    (str): the path to the .onnx file
        vnnlib  (str): the path to the .vnnlib file
        timeout (int): the time (in ms) to wait before proceeding to the next instance
    """
    
    # print("Looking for connections")
    # eng_name = matlab.engine.find_matlab()
    # print(eng_name)
    # try:
    #     eng = matlab.engine.connect_matlab(eng_name[0])
    #     print(eng)
    # except:
    #     print("Connect to anything")
    #     eng = matlab.engine.connect_matlab()
    #     print(eng)

    # print(eng)

    eng = matlab.engine.start_matlab()

    eng.addpath(os.getcwd())
    # NNV root: prefer $NNV_ROOT if it is exported in the environment (an optional override
    # the eval host may set), else derive it relative to this file
    # (code/nnv/examples/Submission/VNN_COMP2026/execute.py -> code/nnv is 3 dirs up).
    # The old hardcoded '/home/ubuntu/toolkit/code/nnv/' broke whenever the eval VM cloned
    # the toolkit anywhere else, silently dropping NNV from the path -> every run errored.
    nnv_root = os.environ.get('NNV_ROOT')
    if not nnv_root:
        here = os.path.dirname(os.path.abspath(__file__))
        nnv_root = os.path.abspath(os.path.join(here, '..', '..', '..'))
    eng.addpath(eng.genpath(nnv_root))
    # CRITICAL: genpath(nnv_root) adds BOTH this submission folder (VNN_COMP2026) AND the frozen
    # sibling VNN_COMP2025. genpath lists them alphabetically (2025 before 2026) and addpath prepends,
    # so VNN_COMP2025/run_vnncomp_instance.m would SHADOW the 2026 one -> the wrong (frozen 2025)
    # harness runs (observed: "VNN_COMP2025/run_vnncomp_instance.m ... ONNX model not supported").
    # os.getcwd() above does NOT fix this (run_instance.sh does not cd into the 2026 dir). addpath
    # prepends, so adding THIS file's own directory LAST puts the 2026 submission at the top of the
    # path and guarantees run_vnncomp_instance + its local helpers resolve to 2026.
    here = os.path.dirname(os.path.abspath(__file__))
    eng.addpath(here)
    print("Paths added (NNV root: %s; prioritized submission dir: %s)" % (nnv_root, here))
    ## eng.addpath(eng.genpath('/root/Documents/MATLAB/SupportPackages/R2024a')) # This is where the support packages get installed from mpm

    # COMPETITION TIMEOUT WIRING: drive NNV's INTERNAL sound caps from the OFFICIAL per-instance timeout
    # (the toolkit reads the benchmark instances.csv and passes it as run_instance.sh $6 -> here), so the
    # Star/PosLin reach and the acas BaB abort to a sound IN-TIME 'unknown' rather than only being hard-
    # killed by future.result(timeout) below (which remains the backstop). These are the SAME env knobs the
    # dev sweeps use -- only the SOURCE of the budget changes (official timeout vs a fixed dev value), so the
    # verification config/algorithm is unchanged. Per-category / NNV-side, NOT per-instance tuning. Set via
    # eng.setenv so the value lands in the MATLAB engine process (os.environ set after start_matlab would not).
    to = float(timeout)
    eng.setenv('NNV_REACH_BUDGET', '%.3f' % max(1.0, 0.95 * to), nargout=0)   # reach/verify cap ~ official timeout
    if 'NNV_ACAS_BAB_TIMECAP' not in os.environ:                              # let an explicit export win
        eng.setenv('NNV_ACAS_BAB_TIMECAP', '%.3f' % min(60.0, max(2.0, 0.5 * to)), nargout=0)  # acas BaB share of the budget

    status = 2 #initialize with an 'Unknown' status
    future = eng.run_vnncomp_instance(category, onnx, vnnlib, outputlocation, nargout = 2, background=True)

    # print("future initiated")

    timeout = float(timeout)
    # print('Trying to get the results without specified timeout')

    # time.sleep(timeout) # wait until timeout everytime?

    # [status, total_time] = future.result()
    
    crashed = False
    try:
        [status, total_time] = future.result(timeout)
        #print('extra time = ',int(toc-tic))
    except matlab.engine.TimeoutError:
        print("timeout")
        #print('extra time = ',int(toc-tic))
        total_time = timeout
        status = 3
    except Exception as e:
        # Any uncaught MATLAB execution error (e.g. a malformed imported net) used to
        # propagate here and leave NO result file -> the harness sees an empty/missing
        # result. A crash means NNV could not decide, so the SOUND verdict is `unknown`
        # (0 points, never a wrong verdict). We never overwrite a verdict the runner
        # already wrote (a late crash after emit), only fill an empty/missing file.
        print("run_vnncomp_instance raised -> sound unknown:", repr(e))
        total_time = timeout
        status = 4
        crashed = True

    future.cancel()
    eng.quit()

    if status == 3:
        resultfile = outputlocation
        with open(resultfile, 'w') as f:
            f.write('timeout')
    elif crashed:
        # write `unknown` ONLY if the runner did not already write a (valid) verdict
        existing = ''
        try:
            with open(outputlocation) as f:
                existing = f.read().strip()
        except Exception:
            existing = ''
        if not existing:
            with open(outputlocation, 'w') as f:
                f.write('unknown')
            print("wrote sound 'unknown' for crashed instance (no result file was produced)")
    # All the other results are written from matlab


def _get_args() -> None:
    """Get the arguments passed to the script from the command line.
    
    Expected usage is : [ACTION, PATH_TO_ONNX, PATH_TO_VNNLIB, TIMEOUT, OUTPUTLOCATION]
    """
    args = sys.argv[1:]
    ACTION = args[0]

    # prepare_instance expects: benchmark_category, onnx, vnnlib
    if (ACTION == 'prepare_instance'):
        if len(args) != 4:
            raise ValueError(f'Incorrect number of arguments, expected 4 got {len(args)}.')
        args.append(None) # timeout
        args.append(None) # outputlocation
    
    # run_instance expects: benchmark_category, onnx, vnnlib, timeout, outputlocation
    if (ACTION == 'run_instance'):
        if len(args) != 6:
            raise ValueError(f'Incorrect number of arguments, expected 6 got {len(args)}.')
    
    print(args)

    return args


if __name__=="__main__":
    # parse the arguments.
    ACTION, CATEGORY, PATH_TO_ONNX, PATH_TO_VNNLIB, TIMEOUT, OUTPUTLOCATION = _get_args()

    # implement logic for each action we might want to take.
    switcher = {
        'prepare_instance': lambda: prepare_instance(CATEGORY, PATH_TO_ONNX, PATH_TO_VNNLIB), # prepare_instance(PATH_TO_ONNX, PATH_TO_VNNLIB),
        'run_instance': lambda: run_instance(CATEGORY, PATH_TO_ONNX, PATH_TO_VNNLIB, TIMEOUT, OUTPUTLOCATION) # run_instance(PATH_TO_ONNX, PATH_TO_VNNLIB, TIMEOUT, OUTPUTLOCATION),
    }

    # retrieve the correct function call based on the input action.
    func = switcher.get(ACTION, 'Invalid')()

    if func == 'Invalid':
        raise ValueError(f'Incorrect ACTION. Expected one of {list(switcher.keys())}; instead got {ACTION}.')