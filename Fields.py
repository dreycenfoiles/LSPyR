from math import pi

from scipy.optimize import minimize
from ngsolve import *

ngsglobals.msg_level = 0


def IncidentWave(wavelength):
	k = 2*pi*n / wavelength
	return CoefficientFunction((0,0,exp(-1J*k*x)))

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
	
	Esc = GetEsc(wavelength)
	Esc_approx = GridFunction(fes)
	Esc_approx.Set(Esc)
	err_func = Norm(Esc-Esc_approx)
	elerr_phy = Integrate(err_func,mesh,element_wise=True,definedon=p)

	return -maxerr

def RefineMesh():

	res = minimize(Error,bounds=(400,800),tol=3)
	wavelength = res.x
	print("Max error wavelength: ",wavelength)
	fesLO = HCurl(mesh,order=1,complex=True)

    for x in range(2):
        print("DoF: ",fes.ndof)
        maxerr = -Error(wavelength)
        print("Max Error: ",maxerr)
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

	ext = 1e-18*k*Integrate((eps_r-1)*E*Conj(Einc),mesh,definedon=mesh.Materials('gold')).imag
	print("Extinction: ",ext)
	return ext





mesh = NewMesh()
fes = HCurl(mesh,order=2,complex=True)
E,W = fes.TnT()

p=pml.BrickRadial((-domain,-domain,-domain),(domain,domain,domain),alpha=1J)
mesh.SetPML(p,'pml')

RefineMesh()

ExtinctionSpectrum(400,700,40)
