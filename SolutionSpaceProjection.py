import scipy.sparse as sp
from scipy.linalg import inv, orth

def SolutionSpaceGeneration(fmin,fmax,wave_list=[]):
	
	global vec_list
	print(wave_list)

	fmid = (fmin + fmax) / 2

	for wavelength in [fmin,fmax]:
		if wavelength not in wave_list:
			Esc = GetEsc(wavelength)
			b = Esc.vec.FV().NumPy()
			vec_list.append(b)
			wave_list.append(wavelength)

	error1 = SolutionSpaceProjection((fmin+fmid)/2,GetError=True)
	error2 = SolutionSpaceProjection((fmid+fmax)/2,GetError=True)
	

	if (error1 < 1e-3) and (error2 < 1e-3):
		return 

	if error1 > error2:
		SolutionSpaceGeneration(fmin,fmid,wave_list=wave_list)
	else:
		SolutionSpaceGeneration(fmid,fmax,wave_list=wave_list)


def SolutionSpaceProjection(wavelength,GetError=False):
	
	global vec_list

	Esc = GridFunction(fes)
	a,f = GetEsc(wavelength,mat=True)
	X = orth(np.column_stack(vec_list))
	Xh = np.conjugate(X.T)
	rows,cols,vals = a.mat.COO()
	A = sp.csr_matrix((vals,(rows,cols)))
	y = f.FV().NumPy()
	z = Xh@y
	v = inv(Xh@A@X)@z

	if GetError == True: 
		error = np.linalg.norm(y - A@X@v)/np.linalg.norm(y)
		return error

	Esc.vec.FV().NumPy()[:] = X@v

	return Esc
