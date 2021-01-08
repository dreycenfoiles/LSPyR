import csv
from math import pi

import matplotlib.pyplot as plt
import numpy as np
from ngsolve import *

from Geometry import *
from Materials import Gold

n = 1.33

ngsglobals.msg_level = 3


class Simulation:

	def __init__(self, geometry, n=1.33):

		self.mesh = geometry[0]
		self.radius = geometry[1]
		self.length = geometry[2]
		self.n = n
		self.fes = HCurl(self.mesh, order=2, complex=True)
		self.fesLO = HCurl(self.mesh, order=1, complex=True)
		self.E, self.W = self.fes.TnT()

	def IncidentWave(self, wavelength):

		k = 2*pi*self.n / wavelength
		return CoefficientFunction((0, exp(-1J*k*x), 0))

	def RelativePermittivity(self, wavelength, delta_n=0.2):

		mesh = self.mesh
		radius = self.radius
		rod_length = self.length

		a = sqrt(x**2+z**2)
		# sqrt(y**2) is for absolute value of y
		r = sqrt(x**2+(sqrt(y**2)-rod_length/2)**2+z**2)
		r0 = sqrt(x**2+y**2+z**2)

		mt_sphere = n + delta_n*(radius/r0)**2
		mt_mid = n + delta_n*(radius/a)
		mt_end = n + delta_n*(radius/r)**2

		permittivities = {"water": 1.33**2,
						  "gold": Gold(wavelength),
						  "pml": 1.33**2,
						  "mt_mid": mt_mid**2,
						  "mt_end": mt_end**2,
						  "mt_sphere": mt_sphere**2,
						  "mt_cyl": 1.414**2}

		return CoefficientFunction([permittivities[mat] for mat in mesh.GetMaterials()])

	def GetEsc(self, wavelength):

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
			mesh.Materials('mt_end') + mesh.Materials('mt_sphere') + \
			mesh.Materials('mt_cyl')

		f = LinearForm(fes)
		f += (eps_r-n**2)*k**2*Einc*W*dx(p)

		c = Preconditioner(a, 'bddc', inverse='pardiso')
		a.Assemble()
		f.Assemble()
		Esc.vec.data = solvers.GMRes(
			a.mat, f.vec, pre=c.mat, printrates=False, restart=750)

		return Esc

	def Error(self, wavelength):

		mesh = self.mesh
		fesLO = self.fesLO

		Esc = self.GetEsc(wavelength)
		Esc_approx = GridFunction(fesLO)
		Esc_approx.Set(Esc)

		p = mesh.Materials('gold') + mesh.Materials('mt_mid') + \
			mesh.Materials('mt_end') + mesh.Materials('mt_sphere') + \
			mesh.Materials('mt_cyl') + mesh.Materials('water')

		err_func = Norm(grad(Esc)-grad(Esc_approx))**2

		# elerr = Integrate(err_func, mesh, element_wise=True)
		elerr = Integrate(err_func, mesh, element_wise=True, definedon=p)
		maxerr = max(elerr)

		return elerr, maxerr

	def RefineMesh(self, fmin=4, fmax=7, tol=1e-1, percentage=.3, iteration=0):

		fes = self.fes
		fesLO = self.fesLO
		mesh = self.mesh

		fmid = (fmin + fmax)/2
		elerr, maxerr = self.Error(fmid)

		print("DoF: ", fes.ndof)
		print("Max Error: ", maxerr)

		p = mesh.Materials('gold') + mesh.Materials('mt_mid') + \
			mesh.Materials('mt_end') + mesh.Materials('mt_sphere') + \
			mesh.Materials('mt_cyl') + mesh.Materials('water')

		if maxerr < tol or iteration > 5:
			return

		else:
			for el in mesh.Elements(p.VB()):
				mesh.SetRefinementFlag(el, elerr[el.nr] > percentage*maxerr)

			mesh.Refine()
			fes.Update()
			fesLO.Update()

			fmid1 = (fmin + fmid)/2
			fmid2 = (fmax + fmid)/2

			error1 = self.Error(fmid1)[1]
			error2 = self.Error(fmid2)[1]

			if error1 > error2:
				self.RefineMesh(fmin=fmid1, fmax=fmid, tol=tol, percentage=percentage, iteration=iteration+1)
			else:
				self.RefineMesh(fmin=fmid, fmax=fmid2, tol=tol, percentage=percentage, iteration=iteration+1)


	def Extinction(self, wavelength):
		print("Wavelength: ", wavelength)

		wavelength /= 100

		mesh = self.mesh

		k = 2*pi / wavelength
		eps_r = self.RelativePermittivity(wavelength)

		Einc = self.IncidentWave(wavelength)
		Esc = self.GetEsc(wavelength)
		E = Einc + Esc

		p = mesh.Materials('gold') + mesh.Materials('mt_mid') + \
			mesh.Materials('mt_end') + mesh.Materials('mt_sphere') + \
			mesh.Materials('mt_cyl')

		ext = 1e-14*k*Integrate((eps_r-1)*E*Conj(Einc), mesh, definedon=p).imag
		return -ext

	def DrawField(self, E):

		vtk = VTKOutput(ma=mesh,
						coefs=[Norm(E)],
						names=["E"],
						filename="result",
						subdivision=2)
		vtk.Do()
