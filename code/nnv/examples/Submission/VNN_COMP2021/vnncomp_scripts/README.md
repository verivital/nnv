1. Tool Name: nnv

2. Manual Setup or Licenses: 

##### Matlab set up with licensing with at least the following toolboxes:

       Computer Vision
       Control Systems
       Deep Learning
       Image Processing
       Optimization
       Parallel Processing
       System Identification
	   
##### installation file

		nnv\code\nnv\dep_requirements.txt
		pip3 install "dep_requirements.txt"
		
##### Add the python matlab engine

		cd /MATLAB/extern/engines/python (this location will be different based on matlab dierctory)
		python3 setup.py install
		
##### Scripts for VNNCOMP2021

	nnv\code\nnv\examples\Submission\VNN_COMP2021\vnncomp_scripts
	
	Scripts are written to run from this location only.
		
3. Use CPU/GPU AWS Instance? CPU

4. Result: 
	
	a. Unknown/Holds: If the given output spec is satisfied,'Unknown/Holds' status depending on reachability method; if 'exact-star' then Holds, if 'approx-star' then Unknown
	
	b. Violates: If the given output spec is not satisfied
	
	c. Timeout: If the code doesn't return any output before mentioned timeout
