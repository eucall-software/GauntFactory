import numpy
import pylab
from scipy import interpolate
from scipy import optimize


power_percent = 100
## Read data from file.
data = numpy.loadtxt('../../daten_zastrau_2015/Bremsspectrum_Laser%dpercent.dat' % (power_percent), delimiter=',')

## Extract data from columns.
wavelength_nm = data[:,0]
dwl = wavelength_nm[1:] - wavelength_nm[:-1]
dwl = dwl.tolist()
dwl.append(dwl[-1])
dwl = numpy.array(dwl)

accepted_range = numpy.where(wavelength_nm < 18.0)
wavelength_nm = wavelength_nm[accepted_range]
accepted_range = numpy.where(wavelength_nm > 8.0)
l_nm = wavelength_nm[accepted_range]
l_nm_fit = numpy.linspace(5,25,201)
d_nm = dwl[accepted_range]
flux = data[:,1]
flux = flux[accepted_range]
norm = sum(flux * d_nm)
normed_flux = flux / norm
length = len(l_nm)

# Determine energy range.
energies_eV = 1239.8/l_nm
energies_Ha = energies_eV/27.2

Ha_to_eV = 27.2
t = numpy.linspace(1,70,70)


spectra_nm = numpy.loadtxt('thermal_bremsstrahlung_T1-70eV.txt').transpose()
n_spectra = spectra_nm.shape[0]
# Multiply by wavelength to convert to photon number spectrum
spectra_nm *= l_nm_fit
d_nm_fit = l_nm_fit[1] -l_nm_fit[0]
# Resample on experimental spectrum:
fit_splines = [interpolate.interp1d(l_nm_fit, spectra_nm[i]) for i in range(n_spectra)]
resampled_spectra = [fit_splines[i](l_nm) for i in range(n_spectra)]

norm = numpy.sum(resampled_spectra*d_nm, axis = 1)
#normalized_spectra = numpy.array([resampled_spectra[i] / norm[i] for i in range(n_spectra)])
normalized_spectra = numpy.array([resampled_spectra[i] for i in range(n_spectra)])

def weights(T1, T2, w1, w2, r):
    w = numpy.exp(-(t-T1)**2/w1) + r*numpy.exp(-(t-T2)**2/w2)
    return w/sum(w)

def weighted_spectra_sum(weights):
    spec = numpy.sum(numpy.array([weights[i] * normalized_spectra[i] for i in range(n_spectra)]), axis=0)
    return spec / sum(spec * d_nm)

def penalty(parameters):
    T1, T2, w1, w2, r = parameters
    if (T1 < min(t) or T1 > max(t)) or \
       (T2 < min(t) or T2 > max(t)) or \
       (w1 < 0. or w2 < 0.) or \
       (r < 0.):
           return 1e6

    return numpy.sum((weighted_spectra_sum(weights(T1, T2, w1, w2, r))- normed_flux)**2)/n_spectra**2

def plot(parameters):
    fitted_spectrum = weighted_spectra_sum(weights(*parameters))
    pylab.plot(l_nm, fitted_spectrum)
    pylab.plot(l_nm, normed_flux)
    #pylab.show()

popt = optimize.minimize(penalty, x0=[10., 30., 1., 1., 1.], method='Nelder-Mead')
bounds=[(10.0, 50.0), (10., 50.), (1.0, 100.), (1., 100.)]
#popt = optimize.minimize(penalty, x0=[10., 30., 1., 1., 1.], method='L-BFGS-B', bounds=bounds)

print popt.x
print penalty(popt.x)

# Plot spectra
plot(popt.x)
pylab.show()

# Plot T distribution.

pylab.plot(t, weights(*popt.x))
pylab.show()
#pylab.plot(l_nm, fitted_spectrum)
#pylab.plot(l_nm, normed_flux)
#pylab.show()
#popt, procc = optimize.minimize( penalty, x0=(10., 30., 1., 1., 1.), method='Nelder-Mead', bounds=[(0.,1.) for i in range(n_spectra)] )




#outdata = spectra_nm
#outdata = numpy.vstack(outdata).transpose()
##print outdata
#numpy.savetxt('thermal_bremsstrahlung_T10-50eV.txt', outdata)
## Normalize to max.
#spectra_nm = [spectrum / spectrum[middle] for spectrum in spectra_nm]

### Plot.
#T = [10,20,30,40,50]
#pylab.semilogy(l_nm, flux, '*')
#[pylab.semilogy(l_nm, spectrum_nm, label='T = %d eV' %(T[i])) for  i,spectrum_nm in enumerate(spectra_nm)]

#pylab.title('%d/100 laser power' % (power_percent))
#pylab.xlabel('Wavelength (nm)')
#pylab.ylabel('Spectral flux')
#pylab.legend(loc=4)
#pylab.show()




