import scipy.sparse as sp
from scipy.linalg import inv, orth
from Simulation import Simulation


class SolutionSpace(Simulation):

	tol = 1e-3

	def __init__(self, fmin, fmax):

		self.Generate(fmin, fmax)
		self.vec_list = []

		for wavelength in [fmin, fmax]:
			Esc = self.GetEsc(wavelength)
			b = Esc.vec.FV().NumPy()
			self.vec_list.append(b)

		fmid = (fmin + fmax) / 2
		error = self.Projection((fmin+fmid)/2, GetError=True)
		if error < tol:
			return
		else:
			Generate(fmin, fmax)

	def Generate(self, fmin, fmax):

		fmid = (fmin + fmax) / 2

		error1 = self.Projection((fmin+fmid)/2, GetError=True)
		error2 = self.Projection((fmid+fmax)/2, GetError=True)

		if (error1 < tol) and (error2 < tol):
			return

		elif error1 >= error2:
			self.Generate(fmin, fmid)
		else:
			self.Generate(fmid, fmax)

	def Projection(self, wavelength, GetError=False):

		Esc = GridFunction(self.fes)
		a, f = self.GetEsc(wavelength, mat=True)
		X = orth(np.column_stack(self.vec_list))
		Xh = np.conjugate(X.T)
		rows, cols, vals = a.mat.COO()
		A = sp.csr_matrix((vals, (rows, cols)))
		y = f.FV().NumPy()
		z = Xh@y
		v = inv(Xh@A@X)@z

		if GetError:
			error = np.linalg.norm(y - A@X@v)/np.linalg.norm(y)
			return error

		else:
			Esc.vec.FV().NumPy()[:] = X@v
			return Esc

	# def Spectrum(self, )

	


