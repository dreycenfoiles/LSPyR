from itertools import cycle
from math import pi

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.interpolate import interp1d
from scipy.linalg import inv, orth
from scipy.optimize import minimize_scalar
from scipy.constants import epsilon_0

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


def Au_Permittivity(wavelength):
	energy = np.array([0.64, 0.77, 0.89, 1.02, 1.14, 1.26, 1.39, 1.51, 1.64, 1.76, 1.88, 2.01, 2.13, 2.26, 2.38, 2.50, 2.63, 2.75, 2.88, 3.00, 3.12, 3.25, 3.37, 3.50, 3.62, 3.74, 3.87, 3.99, 4.12, 4.24, 4.36, 4.49, 4.61, 4.74, 4.86, 4.98, 5.11, 5.23, 5.36, 5.48, 5.60, 5.73, 5.85, 5.98, 6.10, 6.22, 6.35, 6.47, 6.60])

	wavelengths = 1239.84/energy 

	refractiveindex = np.array([0.92 + 13.78j, 0.56 + 11.21j, 0.43 + 9.519j, 
	0.35 + 8.145j, 0.27 + 7.150j, 0.22 + 6.350j, 0.17 + 5.663j, 
	0.16 + 5.083j, 0.14 + 4.542j, 0.13 + 4.103j, 0.14 + 3.697j, 
	0.21 + 3.272j, 0.29 + 2.863j, 0.43 + 2.455j, 0.62 + 2.081j, 
	1.04 + 1.833j, 1.31 + 1.849j, 1.38 + 1.914j, 1.45 + 1.948j, 
	1.46 + 1.958j, 1.47 + 1.952j, 1.46 + 1.933j, 1.48 + 1.895j, 
	1.50 + 1.866j, 1.48 + 1.871j, 1.48 + 1.883j, 1.54 + 1.898j, 
	1.53 + 1.893j, 1.53 + 1.889j, 1.49 + 1.878j, 1.47 + 1.869j, 
	1.43 + 1.847j, 1.38 + 1.803j, 1.35 + 1.749j, 1.33 + 1.688j, 
	1.33 + 1.631j, 1.32 + 1.577j, 1.32 + 1.536j, 1.30 + 1.497j, 
	1.31 + 1.460j, 1.30 + 1.427j, 1.30 + 1.387j, 1.30 + 1.350j, 
	1.30 + 1.304j, 1.33 + 1.277j, 1.33 + 1.251j, 1.34 + 1.226j, 
	1.32 + 1.203j, 1.28 + 1.188j])

	permittivity = refractiveindex**2
	real_permittivity = permittivity.real
	imag_permittivity = permittivity.imag 
	
	real_fit = interp1d(wavelengths,real_permittivity,kind='cubic')
	imag_fit = interp1d(wavelengths,imag_permittivity,kind='cubic')

	return np.complex(real_fit(wavelength),imag_fit(wavelength))

def IncidentWave(wavelength):
	k = 2*pi*n / wavelength
	return CoefficientFunction((0,0,exp(-1J*k*x)))


def RelativePermittivity(wavelength):

	permittivities = {'water' : 1.33**2, 'gold' : Au_Permittivity(wavelength), 'pml' : 1.33**2}

	return CoefficientFunction([permittivities[mat] for mat in mesh.GetMaterials()])


def GetEsc(wavelength):

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

# mesh.ngmesh.Save("nanosphere.vol")

# Esc = GetEsc(550)
# Einc = IncidentWave(550)
# E = Esc + Einc 
# DrawField(E)

ExtinctionSpectrum(400,700,40)
