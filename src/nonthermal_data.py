import numpy
from matplotlib import pylab

from FreeFreeCalculator import FreeFreeCalculator

ni = 1e21 # cm^-3
Z = 1.0

wavelengths = numpy.linspace(7,21,101) #nm
energies = 1239./wavelengths #eV
#print energies

# Read velocity PIC data.
v_data = numpy.loadtxt('../data/Distrib138.dat')
v_axis = v_data[:,0]

dv_raw = v_axis[1:] - v_axis[:-1]
dv = numpy.empty(len(dv_raw)+1)
dv[:-1] = dv_raw
dv[-1] = dv_raw[-1]
p_axis = v_data[:,1] * dv
p_axis = p_axis / sum(p_axis)

# Normalize

ff_calculator = FreeFreeCalculator(ni, Z, energies, v_axis)

spectrum = ff_calculator.freeFreeSpectrum()
spectra = ff_calculator.freeFreeSpectra()

pylab.figure(1)
pylab.plot(energies, spectrum)
pylab.show()

