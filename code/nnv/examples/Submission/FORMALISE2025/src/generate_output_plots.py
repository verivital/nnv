from dotenv import load_dotenv
import matlab.engine
import os

if __name__=="__main__":
    # load environment variables
    load_dotenv()

    # define global variables
    NNV_PATH = os.environ['NNV_PATH']
    NPY_MATLAB_PATH = os.environ['NPY_MATLAB_PATH']

    eng = matlab.engine.start_matlab()
    print('started matlab engine!')

    # add nnv path and npy-matlab path
    eng.addpath(os.getcwd())
    eng.addpath(eng.genpath(NNV_PATH))
    eng.addpath(eng.genpath(NPY_MATLAB_PATH))

    # run the script
    eng.generate_output_plots(nargout=0)

