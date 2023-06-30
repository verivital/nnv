"""Launch MATLAB Engines with Python

This module was created for local testing of reachability analysis (https://github.com/verivital/nnv.git) ahead of VNN-COMP 2023.
To perform testing, the following functionality is enabled.

- Start a MATLAB engine.
- Run an instance of reachability analysis given some arguments through the open MATLAB engine.

Arguments can be passed to the call to the script and will be parsed accordingly.

Created by: Samuel Sasaki (VeriVITAL)
Date: June 28, 2023
"""


import sys
import matlab.engine
import time


def prepare_instance(onnx: str, vnnlib: str) -> None:
    """Set up the MATLAB engine for running an instance.

    Parameters:
        onnx   (str): the path to the .onnx file
        vnnlib (str): the path to the .vnnlib file
    """
    # start matlab engine as a shared engine
    eng = matlab.engine.start_matlab(background=True, option='-r "matlab.engine.shareEngine"')

    # keep MATLAB engine open until manually killed
    while True:
        time.sleep(0.5)
    # eng.quit()


def run_instance(category, onnx, vnnlib, outputlocation) -> None:
    """Run an instance based on parameters defined in .csv file.

    Parameters:
        onnx    (str): the path to the .onnx file
        vnnlib  (str): the path to the .vnnlib file
        timeout (int): the time (in ms) to wait before proceeding to the next instance
    """
    
    # find the matlab engine and connect to it;
    # there should only be 1 running matlab engine.
    eng_name = matlab.engine.find_matlab()[0]
    eng = matlab.engine.connect_matlab(name=eng_name)

    print(f'Successfully connected to engine: {eng_name}.')

    # add logic for running the reachability analysis.

    status = 0 #initialize with an 'Unknown' status
    #toc = time.perf_counter()
    #print('timestep :',toc) 
    future = eng.run_vnncomp_instance(category, onnx, vnnlib, outputlocation,background=True)
    
    try: 
        [status, total_time] = future.result(timeout=args.timeout)
        #print('extra time = ',int(toc-tic))
    except matlab.engine.TimeoutError:
        print("timeout")
        #print('extra time = ',int(toc-tic))
        total_time = args.timeout
        status = 3
        
    future.cancel()
    eng.quit() 
    
    if status == 0:
        result = 'unknown' #Unknown
        #print('Unknown and time: ',total_time)
    elif status == 1:
        result = 'holds'
        #print('Holds and time: ',total_time)
    elif status == 3:
        result = 'timeout'
        #print('Timed Out and time: ',total_time)
    
    resultfile = args.outputlocation
    with open(resultfile, 'w') as f:
        f.write(result)

    eng.quit()


def _get_args() -> None:
    """Get the arguments passed to the script from the command line.
    
    Expected usage is : [ACTION, PATH_TO_ONNX, PATH_TO_VNNLIB, TIMEOUT]
    """
    args = sys.argv[1:]
    ACTION = args[0]

    # prepare_instance expects: benchmark_category, onnx, vnnlib
    if (ACTION == 'prepare_instance'):
        if len(args) != 3:
            raise ValueError(f'Incorrect number of arguments, expected 3 got {len(args)}.')
    
    # run_instance expects: benchmark_category, onnx, vnnlib, timeout
    if (ACTION == 'run_instance'):
        if len(args) != 4:
            raise ValueError(f'Incorrect number of arguments, expected 4 got {len(args)}.')
    
    return args


if __name__=="__main__":
    # parse the arguments.
    ACTION, PATH_TO_ONNX, PATH_TO_VNNLIB, *TIMEOUT = _get_args()
    TIMEOUT = TIMEOUT[0] if (ACTION == 'run_instance') else None

    # implement logic for each action we might want to take.
    switcher = {
        'prepare_instance': lambda: prepare_instance(PATH_TO_ONNX, PATH_TO_VNNLIB), # prepare_instance(PATH_TO_ONNX, PATH_TO_VNNLIB),
        'run_instance': lambda: run_instance(PATH_TO_ONNX, PATH_TO_VNNLIB, TIMEOUT) # run_instance(PATH_TO_ONNX, PATH_TO_VNNLIB, TIMEOUT),
    }

    # retrieve the correct function call based on the input action.
    func = switcher.get(ACTION, 'Invalid')()

    if func == 'Invalid':
        raise ValueError(f'Incorrect ACTION. Expected one of {list(switcher.keys())}; instead got {ACTION}.')