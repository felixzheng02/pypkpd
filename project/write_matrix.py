## Function written to match MATLAB function


import project.all_modules as am

def write_matrix(f, x: am.np.ndarray):
	for i in range(0, x.shape[0]):
		am.fprintf(f,"%6e",x[i, :])
        am.fprintf(f,"\n")
	return
