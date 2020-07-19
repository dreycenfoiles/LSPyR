from math import pi

import matplotlib.pyplot as plt
import numpy as np
from ngsolve import *
from scipy.optimize import minimize_scalar
from Geometry import Nanorod
from Materials import Gold
from SolutionSpaceProjection import SolutionSpace
import csv

n = 1.33

ngsglobals.msg_level = 3


def IncidentWave(wavelength):

	k = 2*pi*n / wavelength 
	return CoefficientFunction((0,exp(-1J*k*x),0))

def RelativePermittivity(wavelength,delta_n=0.2):
	# sqrt(y**2) is for absolute value of y

	a = sqrt(x**2+z**2)
	r = sqrt(x**2+(sqrt(y**2)-rod_length/2)**2+z**2)
	r0 = sqrt(x**2+y**2+z**2)

	mt_sphere = n + delta_n*(radius/r0)**2
	mt_mid = n + delta_n*(radius/a)
	mt_end = n + delta_n*(radius/r)**2

	permittivities = {"water" : 1.33**2, 
					   "gold" : Gold(wavelength), 
					   "pml" : 1.33**2, 
					   "mt_mid" : mt_mid**2, 
					   "mt_end" : mt_end**2, 
					   "mt_sphere" : mt_sphere**2}

	return CoefficientFunction([permittivities[mat] for mat in mesh.GetMaterials()])



def GetEsc(wavelength,mat=False):

	Esc = GridFunction(fes)
	Einc = IncidentWave(wavelength)

	k = 2*pi / wavelength
	eps_r = RelativePermittivity(wavelength)

	a = BilinearForm(fes,symmetric=True)
	a += (curl(E)*curl(W)-k*k*eps_r*E*W)*dx
	a += -1J*k*E.Trace()*W.Trace()*ds('outer')

	p = mesh.Materials('gold') + mesh.Materials('mt_mid') + mesh.Materials('mt_end') + mesh.Materials('mt_sphere')

	f = LinearForm(fes)
	f += (eps_r-n**2)*k**2*Einc*W*dx(p)

	c = Preconditioner(a,'bddc')

	with TaskManager():
		a.Assemble()
		f.Assemble()

	if mat:
		return a,f 
	else:
		Esc.vec.data = solvers.GMRes(a.mat,f.vec,pre=c.mat,printrates=False,restart=750)

	return Esc

def Error(wavelength):
	p = mesh.Materials('gold') + mesh.Materials('mt_mid') + mesh.Materials('mt_end') + mesh.Materials('water') + mesh.Materials('mt_sphere')
	Esc = GetEsc(wavelength)
	Esc_approx = GridFunction(fes)
	Esc_approx.Set(Esc)
	err_func = Norm(Esc-Esc_approx)
	elerr_phy = Integrate(err_func,mesh,element_wise=True,definedon=p)
	maxerr = max(elerr_phy)
	print(maxerr)

	return -maxerr

def RefineMesh():

	res = minimize_scalar(Error,bounds=(500,1500),method='Bounded',tol=8)
	wavelength = res.x
	# print("Max error wavelength: ",wavelength)
	print(fes.ndof)

	Break_Condition = False
	while not Break_Condition: 

		Esc = GetEsc(wavelength)
		Esc_approx = GridFunction(fesLO)
		Esc_approx.Set(Esc)
		err_func = Norm(Esc-Esc_approx)
		elerr = Integrate(err_func,mesh,element_wise=True)
		maxerr = -Error(wavelength)
		for el in mesh.Elements():
			mesh.SetRefinementFlag(el, elerr[el.nr] > .5*maxerr and (el.mat != 'pml'))
		mesh.Refine()
		fes.Update()
		fesLO.Update()

		if maxerr < 1:
			Break_Condition = True

	return wavelength

def Extinction(wavelength):
	k = 2*pi / wavelength
	eps_r = RelativePermittivity(wavelength)
	# print("Wavelength: ",wavelength)

	Einc = IncidentWave(wavelength)
	Esc = GetEsc(wavelength)
	E = Einc + Esc
	
	p = mesh.Materials('gold') + mesh.Materials('mt_mid') + mesh.Materials('mt_end') + mesh.Materials('mt_sphere')

	ext = 1e-18*k*Integrate((eps_r-1)*E*Conj(Einc),mesh,definedon=p).imag
	# print("Extinction: ",ext)
	return -ext

def DrawField(E):

	vtk = VTKOutput(ma=mesh,
				coefs=[Norm(E)],
				names = ["E"],
				filename="result",
				subdivision=2)
	vtk.Do()


def ExtinctionSpectrum(InitWave,FinWave,Steps):

	waverange = np.linspace(InitWave,FinWave,Steps)
	extinction = np.array([Extinction(wavelength) for wavelength in waverange])

	plt.figure()
	plt.plot(waverange,extinction)
	plt.savefig('extinction_rod.png')
	plt.show()


def Saturation(aeff,ratio):

	global mesh 
	global fes 
	global fesLO
	global radius 
	global length 
	global rod_length
	global E,W

	print("aeff: ",aeff)
	print("ratio: ",ratio)

	mt_length_list = np.linspace(0,50,11)
	ext_list = []

	for mt_length in mt_length_list:
		print("length: ",mt_length)

		radius = aeff*(2/(3*ratio-1))**(1/3)
		length = 2 * radius * ratio 
		rod_length = length - radius

		domain = length + 250

		mesh = Nanorod(aeff,ratio,mt_length)
		fes = HCurl(mesh,order=2,complex=True)
		fesLO = HCurl(mesh,order=1,complex=True)
		E,W = fes.TnT()

		p=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
		mesh.SetPML(p,'pml')
	
		maxerr_wavelength = RefineMesh()

		res = minimize_scalar(Extinction,method="Bounded",bounds=(500,1500))
		print("Wavelength: ",res.x)

		ext_list.append(res.x)

	with open("Saturation_Data.csv","a") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow([aeff,ratio,ext_list])
	# return ext_list


aeff_list = np.linspace(10,100,4)
ar_list = np.linspace(1,5,5)

# data = {"aeff="+str(aeff)+",ar="+str(ar) : Saturation(aeff,ar) for aeff in aeff_list for ar in ar_list}

for aeff in aeff_list:
	 for ar in ar_list:
		 Saturation(aeff,ar)

# df = pd.DataFrame(data)
# df.to_excel("Nanorod_Saturation_Data.xlsx")
