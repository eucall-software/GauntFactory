import numpy

from GauntCalculator import freeFreeThermalSpectrum


## Extract data from columns.
l_nm = numpy.linspace(5,25,201)

# Determine energy range.
energies_eV = 1239.8/l_nm
energies_Ha = energies_eV/27.2

Ha_to_eV = 27.2
T = numpy.linspace(10,50,41)
spectra_eV = [freeFreeThermalSpectrum(1.0, t/Ha_to_eV, energies_Ha) for t in T]

### Convert to wavelength spectrum.
spectra_nm = [spectrum_eV / l_nm**2 for spectrum_eV in spectra_eV]
outdata = spectra_nm
outdata = numpy.vstack(outdata).transpose()
##print outdata
numpy.savetxt('thermal_bremsstrahlung_T1-70eV.txt', outdata)
