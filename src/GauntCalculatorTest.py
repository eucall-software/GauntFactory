""" File FreeFreeMatrixTest.py
    Contains unittests for the FreeFreeMatrix object.
"""
import numpy
from matplotlib import pylab

from GauntCalculator import *

alpha = 1./137.

def testGauntFactorFreeFree():

    Z = 1.0
    energies = 10.0**numpy.linspace(-5,8,101) # Ry
    energy = 1.0e5
    velocities = numpy.sqrt(energies*alpha**2) # Ha
    velocity = numpy.sqrt(energy*alpha**2) # Ha

    photon_energies = 10.**numpy.arange(3,7) # Ry
    photon_energy = 1.0e4

    #gff = gauntFactorFF(Z,velocity, photon_energy/2)
    #print gff
    gff = [[gauntFactorFF(Z, v, photon_energy/2.) for v in velocities] for photon_energy in photon_energies]

    [pylab.loglog(energies, gff[i],'k') for i in range(len(photon_energies))]
    pylab.xlabel('Electron initial energy (Ry)')
    pylab.ylabel('free-free Gaunt factor')

    pylab.show()

def testFreeFreeEnsembleSpectrum():
    T = 100.0
    vT = math.sqrt(2.*T)*alpha
    #velocities = vT*10.**numpy.linspace(-1,2, 3)
    velocities = vT*numpy.linspace(0.1, 2, 20)
    weights = MaxwellBoltzmannDistribution(velocities, T)
    energies = T*10.**numpy.linspace(-1,2, 128)
    spectra = numpy.array([weight*numpy.array([freeFreeEmissionSpectrum(1.0, v, e)  for e in energies]) for (weight, v) in zip(weights, velocities)])
    S = numpy.sum(spectra, axis=0)

    spectrum = freeFreeEnsembleSpectrum(1.0, T, energies, velocities, weights)

    [pylab.loglog(energies, spectra[i]) for i in range(len(spectra))]
    pylab.plot(energies, spectrum, '-o')
    pylab.plot(energies, S, '-x')
    #pylab.legend(loc=1)
    pylab.show()

def testFreeFreeEnsembleSpectrumTemperature():
    T = numpy.array([1, 10, 100.0])
    energies = numpy.linspace(1, 10, 10)
    for t in T:

        spectrum = freeFreeThermalSpectrum(1.0, t, energies )

        pylab.semilogy(energies, spectrum, label='T = %.2f' % t)

    pylab.legend(loc=1)
    pylab.show()

def testBoundFreeClassical():
    energies = 10.**numpy.linspace(-4, 3, 101)
    sigma = []
    n = 1
    l = 0
    for en in energies:
        sigma.append(crossSectionBoundFreeClassical(1.0, n, en) )
    pylab.loglog(energies, sigma)
    pylab.show()

def testBoundFreeGauntFactor():
    energies = 10.**numpy.linspace(-4, 3, 101)
    gbf = []
    n = 1
    l = 0
    for en in energies:
        g = abs(gauntFactorBoundFree(1.0, n, l, en))
        print en, g
        gbf.append(g)

    pylab.loglog(energies, gbf)
    pylab.show()
#testGauntFactorFreeFree()
#testFreeFreeEnsembleSpectrum()
#testFreeFreeEnsembleSpectrumTemperature()
testBoundFreeClassical()
