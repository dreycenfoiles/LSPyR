from math import pi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netgen.csg import *
from ngsolve import *
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

n = 1.33

ngsglobals.msg_level = 0


def Nanorod(aeff,ratio,mt_length):

	geo = CSGeometry()

	if mt_length == 0:

		radius = aeff*(2/(3*ratio-1))**(1/3)
		length = 2 * radius * ratio 
		cyl_length = length/2 - radius

		physical_space = length + 100
		domain = physical_space + 150

		sphere1 = Sphere(Pnt(0,-cyl_length,0),radius)
		sphere2 = Sphere(Pnt(0,cyl_length,0),radius)

		cyl1 = Cylinder(Pnt(0,-2*cyl_length,0),Pnt(0,2*cyl_length,0),radius)
		cyl2 = Cylinder(Pnt(0,-2*physical_space,0),Pnt(0,2*physical_space,0),radius+100)
		cyl3 = Cylinder(Pnt(0,-2*domain,0),Pnt(0,2*domain,0),radius+200).bc('outer')

		box2 = OrthoBrick(Pnt(-physical_space,-physical_space,-physical_space),Pnt(physical_space,physical_space,physical_space)) 
		box3 = OrthoBrick(Pnt(-domain,-domain,-domain),Pnt(domain,domain,domain)).bc('outer')

		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))

		middle = cyl1*plane1*plane2

		AuNP = (middle+sphere1+sphere2).mat('gold')
		water = (cyl2*box2 - AuNP).mat('water')
		pmldom = (cyl3*box3 - cyl2*box2).mat('pml')

		geo.Add(AuNP)
		geo.Add(water)
		geo.Add(pmldom)

	else:

		radius = aeff*(2/(3*ratio-1))**(1/3)
		length = 2 * radius * ratio 
		cyl_length = length/2 - radius

		particle = radius + mt_length
		physical_space = length + mt_length + 100
		domain = physical_space + 100

		sphere1 = Sphere(Pnt(0,-cyl_length,0),radius)
		sphere2 = Sphere(Pnt(0,cyl_length,0),radius)
		sphere3 = Sphere(Pnt(0,-cyl_length,0),particle)
		sphere4 = Sphere(Pnt(0,cyl_length,0),particle)

		cyl1 = Cylinder(Pnt(0,-cyl_length,0),Pnt(0,cyl_length,0),radius)
		cyl2 = Cylinder(Pnt(0,-cyl_length,0),Pnt(0,cyl_length,0),particle)
		cyl3 = Cylinder(Pnt(0,-physical_space,0),Pnt(0,physical_space,0),particle+100)
		cyl4 = Cylinder(Pnt(0,-domain,0),Pnt(0,domain,0),particle+200).bc('outer')

		box1 = OrthoBrick(Pnt(-radius,-cyl_length,-radius),Pnt(radius,cyl_length,radius)) 
		box2 = OrthoBrick(Pnt(-particle,-cyl_length,-particle),Pnt(particle,cyl_length,particle)) 
		box3 = OrthoBrick(Pnt(-physical_space,-physical_space,-physical_space),Pnt(physical_space,physical_space,physical_space)) 
		box4 = OrthoBrick(Pnt(-domain,-domain,-domain),Pnt(domain,domain,domain)).bc('outer')

		plane1 = Plane(Pnt(0,-cyl_length,0),Vec(0,-1,0))
		plane2 = Plane(Pnt(0,cyl_length,0),Vec(0,1,0))
		
		middle = cyl1*plane1*plane2
		endcap1 = sphere1-plane1
		endcap2 = sphere2-plane2

		AuNP = (middle+sphere1+sphere2).mat('gold')

		middle_mt = ((cyl2-cyl1)*plane1*plane2).mat('mt_mid')
		endcap1_mt = (sphere3-plane1)-sphere1
		endcap2_mt = (sphere4-plane2)-sphere2
		endcaps_mt = (endcap1_mt + endcap2_mt).mat('mt_end')
		total_body = AuNP + middle_mt + endcaps_mt

		water = (cyl3*box3 - total_body).mat('water')
		pmldom = (cyl4*box4 - cyl3*box3).mat('pml')

		geo.Add(AuNP)
		geo.Add(middle_mt)
		geo.Add(endcaps_mt)
		geo.Add(water)
		geo.Add(pmldom)

	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)
	return mesh


def IncidentWave(wavelength):

	k = 2*pi*n / wavelength 
	return CoefficientFunction((0,0,exp(-1J*k*x)))


def Au_permittivity(wavelength):
	energy = np.array([0.64, 0.77, 0.89, 1.02, 1.14, 1.26, 1.39, 1.51, 1.64, 1.76, 1.88, 2.01, 2.13, 2.26, 2.38, 2.50, 2.63, 2.75, 2.88, 3.00, 3.12, 3.25, 3.37, 3.50, 3.62, 3.74, 3.87, 3.99, 4.12, 4.24, 4.36, 4.49, 4.61, 4.74, 4.86, 4.98, 5.11, 5.23, 5.36, 5.48, 5.60, 5.73, 5.85, 5.98, 6.10, 6.22, 6.35, 6.47, 6.60])

	wavelengths = 1239.84/energy

	refractiveindex = np.array([0.92 + 13.78j, 0.56 + 11.21j, 0.43 + 9.519j, 0.35 + 8.145j, 0.27 + 7.150j, 0.22 + 6.350j, 0.17 + 5.663j, 
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

	permitivity = refractiveindex**2
	real_permitivtiy = permitivity.real
	imag_permitivity = permitivity.imag 

	real_fit = interp1d(wavelengths,real_permitivtiy,kind='cubic')
	imag_fit = interp1d(wavelengths,imag_permitivity,kind='cubic')

	return np.complex(real_fit(wavelength),imag_fit(wavelength))


def RelativePermittivity(wavelength,delta_n=0.2):
	# sqrt(y**2) is for absolute value of y

	a = sqrt(x**2+z**2)
	r = sqrt(x**2+(sqrt(y**2)-rod_length/2)**2+z**2)

	mt_mid = 1.33 + delta_n*(a/radius)
	mt_end = 1.33 + delta_n*(r/radius)**2

	permitivities = {"water" : 1.33**2, "gold" : Au_permittivity(wavelength), "pml" : 1.33**2, "mt_mid" : mt_mid**2, "mt_end" : mt_end**2}

	return CoefficientFunction([permitivities[mat] for mat in mesh.GetMaterials()])



def GetEsc(wavelength):

	Esc = GridFunction(fes)
	Einc = IncidentWave(wavelength)

	k = 2*pi / wavelength
	eps_r = RelativePermittivity(wavelength)

	a = BilinearForm(fes,symmetric=True)
	a += (curl(E)*curl(W)-k*k*eps_r*E*W)*dx
	a += -1J*k*E.Trace()*W.Trace()*ds('outer')

	p = mesh.Materials('gold') + mesh.Materials('mt_mid') + mesh.Materials('mt_end')

	f = LinearForm(fes)
	f += (eps_r-n**2)*k**2*Einc*W*dx(p)

	c = Preconditioner(a,'bddc')

	with TaskManager():
		a.Assemble()
		f.Assemble()

	Esc.vec.data = solvers.GMRes(a.mat,f.vec,pre=c.mat,printrates=False,restart=750)

	return Esc

def Error(wavelength):
	print("Test wavelength: ",wavelength)
	p = mesh.Materials('gold') + mesh.Materials('mt_mid') + mesh.Materials('mt_end') + mesh.Materials('water')
	Esc = GetEsc(wavelength)
	Esc_approx = GridFunction(fes)
	Esc_approx.Set(Esc)
	err_func = Norm(Esc-Esc_approx)
	elerr_phy = Integrate(err_func,mesh,element_wise=True,definedon=p)
	maxerr = max(elerr_phy)

	return -maxerr

def RefineMesh():
	
	res = minimize_scalar(Error,bounds=(500,800),method='Bounded',tol=3)
	wavelength = res.x
	print("Max error wavelength: ",wavelength)

	with TaskManager():
		for x in range(2): 	
			print("DoF:", fes.ndof)
			Esc = GetEsc(wavelength)
			Esc_approx = GridFunction(fes)
			Esc_approx.Set(Esc)
			err_func = Norm(Esc-Esc_approx)
			elerr = Integrate(err_func,mesh,element_wise=True)
			maxerr = -Error(wavelength)
			for el in mesh.Elements():
				mesh.SetRefinementFlag(el, elerr[el.nr] > 0 and (el.mat != 'pml'))
			mesh.Refine()
			fes.Update()
			# print(fes.ndof)
			print("maxerr:",maxerr)
			# if maxerr < 1000:
				# return wavelength

def Extinction(wavelength):
	k = 2*pi / wavelength
	eps_r = RelativePermittivity(wavelength)
	print("Wavelength: ",wavelength)

	Einc = IncidentWave(wavelength)
	Esc = GetEsc(wavelength)
	E = Einc + Esc
	
	p = mesh.Materials('gold') + mesh.Materials('mt_mid') + mesh.Materials('mt_end')

	ext = 1e-18*k*Integrate((eps_r-1)*E*Conj(Einc),mesh,definedon=p).imag
	print("Extinction: ",ext)
	return ext

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
		E,W = fes.TnT()

		p=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
		mesh.SetPML(p,'pml')

		err_wavelength = RefineMesh()

		res = minimize_scalar(lambda wavelength: Extinction(wavelength),method="Bounded",bounds=(400,800),tol=1)
		ext_list.append(res.x)

	return ext_list


# aeff_list = np.linspace(10,100,5)
# ar_list = np.linspace(1.1,5,5)

aeff = 40
ratio = 2
mt_length = 0
mesh = Nanorod(aeff,ratio,0)
radius = aeff*(2/(3*ratio-1))**(1/3)
length = 2 * radius * ratio 
rod_length = length - radius
fes = HCurl(mesh,order=2,complex=True)
E,W = fes.TnT()
RefineMesh()
ExtinctionSpectrum(400,700,40)

# data = {"aeff="+str(aeff)+",ar="+str(ar) : Saturation(aeff,ar) for aeff in aeff_list for ar in ar_list}

# df = pd.DataFrame(data)
# df.to_excel("Nanorod_Saturation_Data.xlsx")
