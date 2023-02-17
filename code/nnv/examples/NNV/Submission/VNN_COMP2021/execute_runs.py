'''
Python script to run matlab commands and record the outputs with a timeout.
author: Neelanjana Pal
'''
import argparse
import matlab.engine
import os
import time


def main():  
    parser = argparse.ArgumentParser(description ='script to run vnncomp2021 instance')
    parser.add_argument('onnxfile')
    parser.add_argument('vnnlibfile')
    parser.add_argument('timeout', type = float, default = 100)
    parser.add_argument('outputlocation')
    parser.add_argument('category')
    args = parser.parse_args()
    
    #tic = time.perf_counter()
   
    # start matlab engine
    eng = matlab.engine.start_matlab()
    eng.addpath(os.getcwd())
    eng.addpath(eng.genpath('../nnv/code'))
   
    status = 0 #initialize with an 'Unknown' status
    #toc = time.perf_counter()
    #print('timestep :',toc) 
    future = eng.run_reachability(args.onnxfile,args.vnnlibfile,args.category,nargout=2,background=True)
    
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
    
if __name__ == "__main__":
    main()
