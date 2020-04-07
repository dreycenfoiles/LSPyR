from Fields import *
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def ExtinctionSpectrum(InitWave,FinWave,Steps):


	waverange = np.linspace(InitWave,FinWave,Steps)

	extinction = np.array([Extinction(wavelength) for wavelength in waverange])

	plt.figure()
	plt.xlabel("Wavelength (nm)")
    plt.ylabel(r'Extinction Cross Section $(m^2)$')
	plt.savefig('extinction_comp_40nm.png')


def DrawField(wavelength):

    Esc = GetEsc(wavelength)

	vtk = VTKOutput(ma=mesh,
				coefs=[Norm(E)],
				names = ["E"],
				filename="Efield",
				subdivision=2)
	vtk.Do()