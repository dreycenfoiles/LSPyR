from itertools import cycle
from math import pi

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize_scalar
from scipy.constants import epsilon_0

from netgen.csg import *
from ngsolve import *

ngsglobals.msg_level = 0

particle = 40
physical_space = 300
domain = 450

n = 1.33

def IncidentWave(wavelength):
	k = 2*pi*n / wavelength
	return CoefficientFunction((0,0,exp(-1J*k*x)))


def RelativePermittivity(wavelength):

	permittivities = {'water' : 1.33**2, 'gold' : Au_Permittivity(wavelength), 'pml' : 1.33**2}

	return CoefficientFunction([permittivities[mat] for mat in mesh.GetMaterials()])


def GetEsc(wavelength):

class Electric_Field(object):

	def __init__(self,mesh,wavelength,refine=True):
		self.mesh = mesh
		self.wavelength = wavelength
		self.refine = refine 
	

	fes = HCurl(self.mesh,order=2,complex=True)
	E,W = fes.TnT()

	brick=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
	self.mesh.SetPML(brick,'pml')

	k = 2*pi / wavelength 

	Einc = CoefficientFunction((0,0,exp(-1J*k*n*x)))
	Esc = GridFunction(fes)
	Einc = IncidentWave(wavelength)

	k = 2*pi / wavelength
	eps_r = RelativePermittivity(wavelength)

	a = BilinearForm(fes,symmetric=True)
	a += (curl(E)*curl(W)-k*k*eps_r*E*W)*dx
	a += -1J*k*E.Trace()*W.Trace()*ds('outer')

	f = LinearForm(fes)
	f += (eps_r-n**2)*k**2*Einc*W*dx(mesh.Materials('gold'))

	c = Preconditioner(a,'bddc')

	with TaskManager():
		a.Assemble()
		f.Assemble()

	Esc.vec.data = solvers.GMRes(a.mat,f.vec,pre=c.mat,printrates=False,restart=750)

	return Esc

def Error(wavelength):
	print("Test wavelength: ",wavelength)
	p = mesh.Materials('gold') + mesh.Materials('water')
	Esc = GetEsc(wavelength)
	Esc_approx = GridFunction(fes)
	Esc_approx.Set(Esc)
	err_func = Norm(Esc-Esc_approx)
	elerr_phy = Integrate(err_func,mesh,element_wise=True,definedon=p)
	maxerr = max(elerr_phy)

	return -maxerr

def RefineMesh():

	res = minimize_scalar(Error,bounds=(400,800),method='Bounded',tol=3)
	wavelength = res.x
	print("Max error wavelength: ",wavelength)
	fesLO = HCurl(mesh,order=1,complex=True)

	with TaskManager():
		while True: 	
			p = mesh.Materials('gold') + mesh.Materials('water')
			print("DoF: ",fes.ndof)
			maxerr = -Error(wavelength)
			print("Max Error: ",maxerr)
			p = mesh.Materials('gold') + mesh.Materials('water')
			Esc = GetEsc(wavelength)
			Esc_approx = GridFunction(fesLO)
			Esc_approx.Set(Esc)
			err_func = Norm(Esc-Esc_approx)
			elerr = Integrate(err_func,mesh,element_wise=True)
			for el in mesh.Elements():
				mesh.SetRefinementFlag(el, elerr[el.nr] > .5*maxerr and (el.mat != 'pml'))
			mesh.Refine()
			fes.Update()
			fesLO.Update()


def Extinction(wavelength):

	k = 2*pi / wavelength
	eps_r = RelativePermittivity(wavelength)
	print("Wavelength: ",wavelength)

	Einc = IncidentWave(wavelength)
	Esc = GetEsc(wavelength)
	E = Einc + Esc

	ext = 1e-18*k*Integrate((eps_r-1)*E*Conj(Einc),mesh,definedon=mesh.Materials('gold'),order=10).imag
	print("Extinction: ",ext)
	return ext


def DrawField(E):

	vtk = VTKOutput(ma=mesh,
				coefs=[Norm(E)],
				names = ["E"],
				filename="Efield",
				subdivision=2)
	vtk.Do()

def ExtinctionSpectrum(InitWave,FinWave,Steps):

	from compare2 import data

	print(fes.ndof)

	waverange = np.linspace(InitWave,FinWave,Steps)

	extinction = np.array([Extinction(wavelength) for wavelength in waverange])

	plt.figure()
	plt.plot([el[0]*1e+9 for el in data],[el[1] for el in data],label="Mie")
	plt.plot(waverange,extinction,label="FEM")
	plt.legend()
	plt.savefig('extinction_comp_40nm.png')



mesh = NewMesh()
fes = HCurl(mesh,order=2,complex=True)
E,W = fes.TnT()

p=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
mesh.SetPML(p,'pml')

RefineMesh()

ExtinctionSpectrum(400,700,40)
