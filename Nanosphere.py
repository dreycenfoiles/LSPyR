from itertools import cycle
from math import pi

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.linalg import inv, orth
from scipy.optimize import minimize_scalar
from Materials import *

from netgen.csg import *
from ngsolve import *

ngsglobals.msg_level = 0

particle = 40
physical_space = 300
domain = 450

n = 1.33


def NewMesh():

	geo = CSGeometry()

	sphere1 = Sphere(Pnt(0,0,0),particle)
	sphere2 = Sphere(Pnt(0,0,0),physical_space)
	sphere3 = Sphere(Pnt(0,0,0),domain).bc('outer')

	AuNP = sphere1.mat('gold')
	water = (sphere2 - sphere1).mat('water')
	pml = (sphere3 - sphere2).mat('pml').maxh(80)

	geo.Add(AuNP)
	geo.Add(water)
	geo.Add(pml)

	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)

	return mesh

class Electric_Field(object):

	def __init__(self,mesh,wavelength,refine=True):
		self.mesh = mesh
		self.wavelength = wavelength
		self.Refine = refine 
	
	fes = HCurl(mesh,order=2,complex=True)
	E,W = fes.TnT()

	p=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
	mesh.SetPML(p,'pml')

	k = 2*pi / wavelength 
	Einc = CoefficientFunction((0,0,exp(-1J*k*x)))

	def RelativePermittivity(self,wavelength):

		permittivities = {'water' : 1.33**2, 'gold' : Gold(wavelength), 'pml' : 1.33**2, 'air' : 1}

		return CoefficientFunction([permittivities[mat] for mat in mesh.GetMaterials()])

	def Error(self,wavelength):

		print("Test wavelength: ",wavelength)
		p = mesh.Materials('gold') + mesh.Materials('water')
		Esc = GetEsc(wavelength)
		Esc_approx = GridFunction(fes)
		Esc_approx.Set(Esc)
		err_func = Norm(Esc-Esc_approx)
		elerr_phy = Integrate(err_func,mesh,element_wise=True,definedon=p)
		maxerr = max(elerr_phy)

		return -maxerr

	def RefineMesh(self):
		
		res = minimize_scalar(Error,bounds=(400,700),method='Bounded',tol=3)
		wavelength = res.x
		print("Max error wavelength: ",wavelength)

		with TaskManager():
			while True: 	
				p = mesh.Materials('gold') + mesh.Materials('water')
				print("DoF: ",fes.ndof)
				maxerr = -Error(wavelength)
				print("Max Error: ",maxerr)
				p = mesh.Materials('gold') + mesh.Materials('water')
				Esc = GetEsc(wavelength)
				Esc_approx = GridFunction(fes)
				Esc_approx.Set(Esc)
				err_func = Norm(Esc-Esc_approx)
				elerr = Integrate(err_func,mesh,element_wise=True)
				for el in mesh.Elements():
					mesh.SetRefinementFlag(el, elerr[el.nr] > .5*maxerr and (el.mat != 'pml'))
				mesh.Refine()
				fes.Update()
				if maxerr < 1000:
					return 

	def Calculate(self):

		Esc = GridFunction(fes)
		
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

	def Extinction(self):

		print("Wavelength: ",wavelength)
		
		ext = 1e-18*k*Integrate((eps_r-1)*E*Conj(Einc),mesh,definedon=mesh.Materials('gold')).imag
		print("Extinction: ",ext)
		return ext

	def DrawField(self,E):

		vtk = VTKOutput(ma=mesh,
					coefs=[Norm(E)],
					names = ["E"],
					filename="Efield",
					subdivision=2)
		vtk.Do()

			

def ExtinctionSpectrum(InitWave,FinWave,Steps,name="Extinction.png"):


	waverange = np.linspace(InitWave,FinWave,Steps)

	extinction = np.array([Extinction(wavelength) for wavelength in waverange])

	plt.figure()
	plt.plot(waverange,extinction)
	plt.legend()
	plt.savefig(name)

