# ==============================
# This functions converts the 
# Zonotope approximations of the 
# reach sets and saves them as 
# a struct to use in MATLAB
# ==============================

using MAT

function convert_sets(R, out_file)
	# Get size of generators from set
	setEx = size(Array(R.F.Xk[2].X.generators))
	matGen = zeros(length(R), setEx[1], setEx[2])
	matCenter = zeros(length(R), setEx[1])
	timeZ = zeros(length(R))
	for ss in 1:length(R)
		matGen[ss,:,:] = Array(R.F.Xk[ss].X.generators)
		matCenter[ss,:] = Array(R.F.Xk[ss].X.center)
		timeZ[ss] = R.F.Xk[ss].Δt.hi
	end
	# matwrite(out_file, Dict("matGen" => matGen,"matCenter" => matCenter, "timeV" => timeZ); compress = false)
	matwrite(out_file, Dict("gens" => matGen, "centers" => matCenter, "timeZ" => timeZ); compress = false)
# Dict("gens" => a1, "mats" => a2, "time" => a3)
	return 
end

function convert_boxs(R, out_file)
	matGen = zeros(length(R), 2)
	matCenter = zeros(length(R), 2)
	timeZ = zeros(length(R))
	for ss in 1:length(R)
		matGen[ss,:] = Array(R.F.Xk[ss].X.radius)
		matCenter[ss,:] = Array(R.F.Xk[ss].X.center)
		timeZ[ss] = R.F.Xk[ss].Δt.hi
	end
	matwrite(out_file, Dict("gens" => matGen, "centers" => matCenter, "timeZ" => timeZ); compress = false)
	return
end