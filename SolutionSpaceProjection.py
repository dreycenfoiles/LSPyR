import scipy.sparse as sp
from scipy.linalg import orth, inv
import numpy as np
from Simulation import *


class SolutionSpace(Simulation):

	def __init__(self, mesh, fmin, fmax, refine_tol=1000, percentage=.5, ssp_tol=1e-3, vec_list=[]):
		Simulation.__init__(self, mesh)
		self.RefineMesh(fmin, fmax, tol=refine_tol, percentage=percentage)

		self.vec_list = vec_list
		self.ssp_tol = ssp_tol

		for wavelength in [fmin, fmax]:
			self.Add(wavelength)

		fmid = (fmin + fmax) / 2
		error = self.Projection(fmid, GetError=True)
		print("Initial Error:", error)
		if error < ssp_tol:
			return
		else:
			self.Generate(fmin, fmax)

	def Add(self, wavelength):
		Esc = self.GetEsc(wavelength)
		b = Esc.vec.FV().NumPy()
		self.vec_list.append(b)

	def Generate(self, fmin, fmax):

		fmid = (fmin + fmax) / 2

		error1 = self.Projection((fmin+fmid)/2, GetError=True)
		error2 = self.Projection((fmid+fmax)/2, GetError=True)

		print("Error 1:", error1)
		print("Error 2:", error2)

		if (error1 < self.ssp_tol) and (error2 < self.ssp_tol):
			return

		else:
			self.Add(fmid)

			if error1 >= error2:
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
		y = f.vec.FV().NumPy()
		z = Xh@y
		v = inv(Xh@A@X)@z

		if GetError:
			error = np.linalg.norm(y - A@X@v)/np.linalg.norm(y)
			return error

		else:
			Esc.vec.FV().NumPy()[:] = X@v
			return Esc
	
	def Extinction(self, wavelength):

		mesh = self.mesh

		k = 2*pi / wavelength
		eps_r = self.RelativePermittivity(wavelength)

		Einc = self.IncidentWave(wavelength)
		Esc = self.Projection(wavelength)
		E = Einc + Esc

		p = mesh.Materials('gold') + mesh.Materials('mt_mid') + \
			mesh.Materials('mt_end') + mesh.Materials('mt_sphere')

		ext = 1e-18*k*Integrate((eps_r-1)*E*Conj(Einc), mesh, definedon=p).imag
		return -ext