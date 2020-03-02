#%%
from ngsolve import *
from netgen.csg import *
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
#%%

ngsglobals.msg_level = 2

def Nanorod(length,ratio):

	geo = CSGeometry()

	radius = (length / ratio) / 2
	cyl_end = length/2 - radius 

	cyl1 = Cylinder(Pnt(0,-length,0),Pnt(0,length,0),radius)
	box1 = OrthoBrick(Pnt(-radius,-length,-radius),Pnt(radius,length,radius)) 
	sphere1 = Sphere(Pnt(0,-length,0),radius)
	sphere2 = Sphere(Pnt(0,length,0),radius)
	sphere3 = Sphere(Pnt(0,0,0),15)
	sphere4 = Sphere(Pnt(0,0,0),40).bc('outer')


	AuRod = (cyl1*box1).mat('gold').maxh(1)
	endcap1 = (sphere1 - box1).mat('gold').maxh(1)
	endcap2 = (sphere2 - box1).mat('gold').maxh(1)
	rod = AuRod*endcap1*endcap2
	water = (sphere3 - rod).mat('water').maxh(4)
	pml = (sphere4 - sphere3).mat('pml').maxh(3)


	geo.Add(AuRod)
	geo.Add(endcap1)
	geo.Add(endcap2)
	geo.Add(water)
	geo.Add(pml)

	ngmesh = geo.GenerateMesh()

	mesh = Mesh(ngmesh)
	mesh.Curve(5)
	return mesh


def Au_permitivity(wavelength):
	energy = np.array([0.64, 0.77, 0.89, 1.02, 1.14, 1.26, 1.39, 1.51, 1.64, 1.76, 1.88, 2.01, 2.13, 2.26, 2.38, 2.50, 2.63, 2.75, 2.88, 3.00, 3.12, 3.25, 3.37, 3.50, 3.62, 3.74, 3.87, 3.99, 4.12, 4.24, 4.36, 4.49, 4.61, 4.74, 4.86, 4.98, 5.11, 5.23, 5.36, 5.48, 5.60, 5.73, 5.85, 5.98, 6.10, 6.22, 6.35, 6.47, 6.60])

	wavelengths = 123.984/energy

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

	real_fit = interp1d(wavelengths,real_permitivtiy,kind='linear')
	imag_fit = interp1d(wavelengths,imag_permitivity,kind='linear')

	return np.complex(real_fit(wavelength),imag_fit(wavelength))


def RelativePermitivity(wavelength):

	permitivities = {"water" : 1.33**2, "gold" : Au_permitivity(wavelength), "pml" : 1.33**2}

	return CoefficientFunction([permitivities[mat] for mat in mesh.GetMaterials()])


def GetEfield(wavelength):

	fes.Update()

	Efield = GridFunction(fes)

	k = 2*np.pi/wavelength

	Einc = CoefficientFunction((0,exp(1J*k*x),0))
	curl_Einc = CoefficientFunction((0,0,1J*k*exp(1J*k*x)))

	eps_r = RelativePermitivity(wavelength)

	U = Cross(n_hat,curl_Einc) + 1J*k*Cross(n_hat,Cross(n_hat,Einc))

	a = BilinearForm(fes,symmetric=True)

	a += (curl(E)*curl(Etest)-k*k*eps_r*E*Etest)*dx 
	a += -1J*k*Cross(n_hat,E.Trace())*Cross(n_hat,Etest.Trace())*ds('outer')


	f = LinearForm(fes)
	f += -Etest.Trace()*(Cross(n_hat,curl_Einc) + 1J*k*Cross(n_hat,Cross(n_hat,Einc)))*ds('outer')

	c = Preconditioner(a, 'bddc')
	a.Assemble()
	f.Assemble()
	c.Update()   
	 
	inv = GMRESSolver(mat=a.mat, pre=c.mat, printrates=False) # Use Inverse as preconditioner for next wavelength
	solution = inv * f.vec
	Efield.vec.data = solution
	return Efield


def GetExtinction(wavelength):

	eps_0 = 8.854e-36

	k = 2*np.pi / wavelength

	E = GetEfield(wavelength)

	eps_r = RelativePermitivity(wavelength)
	
	Einc = CoefficientFunction((0,exp(1J*k*x),0))

	P = (eps_r-1)*eps_0*E 

	return 1e+16*4*np.pi*k*Integrate((Conj(Einc)*P).imag,mesh,definedon=mesh.Materials('gold'))


def DrawField(E):

	vtk = VTKOutput(ma=mesh,
				coefs=[Norm(E)],
				names = ["E"],
				filename="result",
				subdivision=2)
	vtk.Do()


def ExtinctionSpectrum(InitWave,FinWave,Steps):

	waverange = np.linspace(InitWave/10,FinWave/10,Steps)
	extinction = np.array([GetExtinction(wavelength) for wavelength in waverange])

	plt.figure()
	plt.plot(waverange*10,extinction)
	plt.savefig('extinction.png')
	plt.show()

	

#%%


mesh = Nanorod(4,1.6)

fes = HCurl(mesh,order=3,complex=True)

E= fes.TrialFunction()
Etest = fes.TestFunction()

p=pml.BrickRadial((-40,-40,-40),(40,40,40),alpha=3J)
mesh.SetPML(p,definedon=mesh.Materials('pml'))

n = specialcf.normal(3)
n_mag = Norm(n)
n_hat = CoefficientFunction((n[0]/n_mag,n[1]/n_mag,n[2]/n_mag))


#%%

# res = minimize_scalar(GetExtinction,bounds=(40,55),method='bounded')
# print(res.x)


#%%
# ExtinctionSpectrum(400,700,15)
E = GetEfield(60)
Draw(Norm(E),mesh,'E')
