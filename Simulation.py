import csv
from math import pi

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from ngsolve import *
from scipy.linalg import inv, orth
from scipy.optimize import minimize_scalar

from Geometry import *
from Materials import Gold

n = 1.33

ngsglobals.msg_level = 3


class Simulation:

	def __init__(self, mesh, n=1.33):

		self.mesh = mesh
		self.n = n
		self.fes = HCurl(mesh, order=2, complex=True)
		self.fesLO = HCurl(mesh, order=1, complex=True)
		self.E, self.W = self.fes.TnT()

	def IncidentWave(self, wavelength):

		k = 2*pi*self.n / wavelength
		return CoefficientFunction((0, exp(-1J*k*x), 0))

	def RelativePermittivity(self, wavelength, delta_n=0.2):

		mesh = self.mesh
		# sqrt(y**2) is for absolute value of y

		# a = sqrt(x**2+z**2)
		# r = sqrt(x**2+(sqrt(y**2)-rod_length/2)**2+z**2)
		# r0 = sqrt(x**2+y**2+z**2)

		# mt_sphere = n + delta_n*(radius/r0)**2
		# mt_mid = n + delta_n*(radius/a)
		# mt_end = n + delta_n*(radius/r)**2

		# permittivities = {"water" : 1.33**2,
		# 				   "gold" : Gold(wavelength),
		# 				   "pml" : 1.33**2,
		# 				   "mt_mid" : mt_mid**2,
		# 				   "mt_end" : mt_end**2,
		# 				   "mt_sphere" : mt_sphere**2,
		# 				   "mt_cyl" : 1.414**2}

		permittivities = {"water" : 1.33**2,
						   "gold" : Gold(wavelength),
						   "pml" : 1.33**2}
		

		return CoefficientFunction([permittivities[mat] for mat in mesh.GetMaterials()])

	def GetEsc(self, wavelength, mat=False):

		mesh = self.mesh
		fes = self.fes
		E = self.E
		W = self.W

		Esc = GridFunction(fes)
		Einc = self.IncidentWave(wavelength)

		k = 2*pi / wavelength
		eps_r = self.RelativePermittivity(wavelength)

		a = BilinearForm(fes, symmetric=True)
		a += (curl(E)*curl(W)-k*k*eps_r*E*W)*dx
		a += -1J*k*E.Trace()*W.Trace()*ds('outer')

		p = mesh.Materials('gold') + mesh.Materials('mt_mid') + \
			mesh.Materials('mt_end') + mesh.Materials('mt_sphere')

		f = LinearForm(fes)
		f += (eps_r-n**2)*k**2*Einc*W*dx(p)

		if mat:
			a.Assemble()
			f.Assemble()
			return a, f
		else:
			c = Preconditioner(a, 'bddc', inverse="sparsecholesky")
			a.Assemble()
			f.Assemble()
			Esc.vec.data = solvers.GMRes(
				a.mat, f.vec, pre=c.mat, printrates=False, restart=750)

			return Esc

	def Error(self, wavelength):

		mesh = self.mesh
		fes = self.fes
		fesLO = self.fesLO

		p = mesh.Materials('gold') + mesh.Materials('mt_mid') + mesh.Materials(
			'mt_end') + mesh.Materials('water') + mesh.Materials('mt_sphere')

		Esc = self.GetEsc(wavelength)
		Esc_approx = GridFunction(fesLO)
		Esc_approx.Set(Esc)

		err_func = Norm(Esc-Esc_approx)
		elerr_phy = Integrate(err_func, mesh, element_wise=True, definedon=p)
		elerr = Integrate(err_func, mesh, element_wise=True)
		maxerr = max(elerr_phy)

		return elerr, maxerr

	def RefineMesh(self, fmin, fmax, tol=1000, percentage=.5):

		fes = self.fes
		fesLO = self.fesLO
		mesh = self.mesh

		fmid = (fmin + fmax)/2
		elerr, maxerr = self.Error(fmid)

		print("DoF: ", fes.ndof)
		print("Max Error: ", maxerr)

		if maxerr < tol:
			return

		else:
			for el in mesh.Elements():
				mesh.SetRefinementFlag(
					el, elerr[el.nr] > percentage*maxerr and (el.mat != 'pml'))

			mesh.Refine()
			fes.Update()
			fesLO.Update()

			fmid1 = (fmin + fmid)/2
			fmid2 = (fmax + fmid)/2

			error1 = self.Error(fmid1)[1]
			error2 = self.Error(fmid2)[1]

			if error1 > error2:
				self.RefineMesh(fmid1, fmid)
			else:
				self.RefineMesh(fmid, fmid2)

	def Extinction(self, wavelength):

		mesh = self.mesh

		k = 2*pi / wavelength
		eps_r = self.RelativePermittivity(wavelength)

		Einc = self.IncidentWave(wavelength)
		Esc = self.GetEsc(wavelength)
		E = Einc + Esc

		p = mesh.Materials('gold') + mesh.Materials('mt_mid') + \
			mesh.Materials('mt_end') + mesh.Materials('mt_sphere')

		ext = 1e-18*k*Integrate((eps_r-1)*E*Conj(Einc), mesh, definedon=p).imag
		return -ext

	def DrawField(self, E):

		vtk = VTKOutput(ma=mesh,
						coefs=[Norm(E)],
						names=["E"],
						filename="result",
						subdivision=2)
		vtk.Do()

	# def ExtinctionSpectrum(self, InitWave, FinWave, Steps):

	#     waverange = np.linspace(InitWave, FinWave, Steps)
	#     extinction = np.array([Extinction(wavelength)
	#                            for wavelength in waverange])

	#     plt.figure()
	#     plt.plot(waverange, extinction)
	#     plt.show()

	# def Saturation(self, aeff, ratio):

	#     print("aeff: ", aeff)
	#     print("ratio: ", ratio)

	#     mt_length_list = np.linspace(0, 50, 11)
	#     ext_list = []

	#     for mt_length in mt_length_list:
	#         print("length: ", mt_length)

	#         radius = aeff*(2/(3*ratio-1))**(1/3)
	#         length = 2 * radius * ratio
	#         rod_length = length - radius

	#         domain = length + 250

	#         mesh = Nanorod(aeff, ratio, mt_length)
	#         fes = HCurl(mesh, order=2, complex=True)
	#         fesLO = HCurl(mesh, order=1, complex=True)
	#         E, W = fes.TnT()

	#         p = pml.BrickRadial((-domain, -domain, -domain),
	#                             (domain, domain, domain), alpha=1J)
	#         mesh.SetPML(p, 'pml')

	#         maxerr_wavelength = RefineMesh()
	#         # SSP = SolutionSpace(maxerr_wavelength-100,maxerr_wavelength+100)

	#         res = minimize_scalar(Extinction, method="Bounded", bounds=(
	#             maxerr_wavelength-100, maxerr_wavelength+100))
	#         print("Wavelength: ", res.x)

	#         ext_list.append(res.x)

	#     with open("Saturation_Data.csv", "a") as csvfile:
	#         writer = csv.writer(csvfile)
	#         writer.writerow([aeff, ratio, ext_list])
		# return ext_list


