""" File FreeFreeMatrixTest.py
    Contains unittests for the FreeFreeMatrix object.
"""
from FreeFreeMatrix import FreeFreeCalculator
import numpy

ni = 1e21
Z = 1
energies = numpy.linspace(1,100,1)/27.2 # 1-100 eV in units of Ha
velocities = [0.1]
ff_calculator = FreeFreeCalculator(ni, Z, energies, velocities)

#self.assertIsInstance(ff_calculator, FreeFreeCalculator)

spectrum = ff_calculator.freeFreeSpectrum()

print spectrum

