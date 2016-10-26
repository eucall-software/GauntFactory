""" File that hosts the free-free coupling matrix L, see Brussaard1962, eq. (3) """
import numpy
import sys
import math
from scipy.constants import physical_constants as PC
from scipy.constants import alpha
import scipy.integrate
from mpmath import hyp2f1
from mpmath import factorial
from mpmath import gamma

Ha_in_eV =27.21138505
aB = PC['Bohr radius']
aB_in_cm = PC['Bohr radius'][0]*100.0

class GauntCalculator:
    def __init__(self,
                 ion_density,
                 ion_charge,
                 energies,
                 velocities,
                 weights=None,
                 ):
        """ Constructor for the FreeFreeCalculator.

        @param ion_density : The number density of ions in units of cm^-3
        @param ion_charge  : The charge state of the ions.
        @param energies    : The energies for which to calculate the free-free spectrum in units of eV.
        @param velocities  : The velocities for which to calculate the radiation in units of c.
        @param weights     : The weights of the velocity distribution.

        """

        self.__ni_cm3 = ion_density
        self.__ni_aB3 = ion_density*(aB_in_cm)**3
        self.__Z = ion_charge

        self.__energies_eV = energies
        self.__energies_Ha = energies/Ha_in_eV
        self.__velocities = velocities

        max_velocity = numpy.max(velocities)
        max_kin_energy = 0.5*max_velocity**2 / alpha**2 * Ha_in_eV
        print 'Maximum kinetic energy = ', max_kin_energy, ' eV'

        min_photon_energy = numpy.min(self.__energies_eV)
        print 'Minimum photon energy in spectral range = ', min_photon_energy, ' eV'

        if max_kin_energy < min_photon_energy:
            print 'WARNING: No free-free emission kinematically allowed.'
            sys.exit()

        de_matrix = numpy.empty((len(self.__velocities), len(self.__energies_Ha)) )
        self.__spectra_v = numpy.empty((len(self.__velocities), len(self.__energies_Ha)), dtype=float )
        for i,v in enumerate(self.__velocities):
            for j,en in enumerate(self.__energies_Ha):
                de_matrix[i,j] = 0.5*v**2/alpha**2 - en
                self.__spectra_v[i,j] = self.__ni_aB3**2 * self.__Z * calculateFreeFreeSpectrum(self.__Z, en, v).real
        if weights is None:
            weights = numpy.ones_like(self.__velocities)
        self.__weights = weights

    def freeFreeSpectrum(self):
        """ Query for the free-free spectrum. """
	return self.__weights.transpose().dot(self.__spectra_v) / sum(self.__weights)

    def freeFreeSpectra(self, weighted=False):
        """ Query for the free-free spectra (for all velocities). """
	if weighted:
	    return [self.__weights[i] * self.__spectra_v[i] / sum(self.__weights) for i in range(len(self.__weights)) ]
	else:
	    return self.__spectra_v

def freeFreeThermalSpectrum(Z, T, energies):
    """  Calculate the emission from an ensemble of electrons

    @param Z : The ion charge.
    @param T : The temperature (Hartree)
    @param enenergies : The energies at which to calculate the spectrum.
    """
    def distribution(e):
        return math.sqrt(e) * math.exp(-e/T)

    norm = scipy.integrate.quad(distribution, 0.0, numpy.inf)[0]
    spectrum = numpy.zeros_like(energies)
    for i,Ephoton in enumerate(energies):
        def integrand(Ekin):
            return distribution(Ekin) * freeFreeEmissionSpectrum(Z, Ephoton, Ekin)
        avg = scipy.integrate.quad(integrand, Ephoton, numpy.inf)[0]
        spectrum[i] =  avg/norm
    return spectrum

def freeFreeEmissionSpectrum(Z, en, Ev):
    """ Calculate the free-free spectrum for a given energy and velocity.
    @param en : Energy in the spectrum for which to calculate the emission.
    @unit     : Ha
    @param Ev : Kinetic energy
    @unit     : Ha
    """

    #print eta_f
    if en >= Ev:
        return 0
    v_initial = alpha * math.sqrt(2.*(Ev - en))
    g_ff = gauntFactorFF(Z, v_initial, en)
    #
    return 64.*math.pi**2/3./math.sqrt(3.) * Z**2 * alpha**3 / math.sqrt(2.*Ev) * g_ff # [Ry/s/Hz*aB^3]
    # cf. Brussaard1962, rearranged into a form found in Karzas1961, eq. 23. note a factor 2 difference in the latter.
    # To get emission spectrum in W/eV/cm^3, multiply by Ry in J , divide by h*1Hz in eV, and divide by (aB in cm)^3.
def G(eta1, eta2, k1, k2, l):
    Gl = abs((k2 - k1) / (k2 + k1))**(1j*eta1 + 1j*eta2) * hyp2f1(l+1-1j*eta2, l+1-1j*eta1, 2.*l+2, -4.*k1*k2/(k1-k2)**2)

    return Gl.real


def I(eta1, eta2, k1, k2, l):

    Il =  0.25*(4.* k1 * k2 / (k1 - k2)**2)**(l+1) * math.exp(math.pi*abs(eta1 - eta2)/2.) * abs(gamma(l+1+1j*eta1) * gamma(l+1+1j*eta2) ) / gamma(2.*l + 2)

    return Il * G(eta1, eta2, k1, k2, l)

def gauntFactorFF_BvdH(Z, electron_velocity, photon_energy):
    """ Free-free Gaunt factor according to Brussaard-van de Hulst (1962), eq. 2. """

    v = electron_velocity
    en = photon_energy

    #eta = Z/sqrt(E/Ry)  = Z/sqrt(0.5mv**2/Ry) = Z/sqrt(mv**2/Ha) = Z/c/sqrt(v**2/c**2) = Z*alpha/v
    eta1 = Z * alpha / v

    #E = hbar^2 k ^2 / 2m =(Ha) (k aB)^2/2
    #=> k = sqrt(2 E) = sqrt( 2* 0.5*mv^2 )  = sqrt(v^2/c^2)*c = v/alpha

    k1 = v / alpha
    # 0.5*mvf**2 = 0.5*mvi**2 - en
    #mvf^2 = mvi^2 + 2*en
    #vf^2 = vi^2 + 2en/m
    #vf = sqrt(vi^2 + 2en/m)
    #vf/c = sqrt(vi^2/c^2 + 2*en/c^2) = sqrt(vi^2 + 2*en*alpha^2)
    eta2 = Z * alpha / numpy.sqrt(v**2 + 2.*en * alpha**2)


    return math.sqrt(3.)/math.pi * L( eta2, eta1)

def L(eta1, eta2):
    """ The L function in Brussaard-van de Hulst (1962), eq. (4). """

    tmp1 = math.pi**2 * eta1 * eta2 / ((math.exp(2.*math.pi*eta1) - 1.) * (1.-math.exp(-2.*math.pi*eta2)))
    tmp2 = abs(DeltaHypergeometric(1j*eta1, 1j*eta2)) / (eta2 - eta1)

    return tmp1 * tmp2

def gauntFactorFF(Z, electron_velocity, photon_energy, use_BvdH=False):
    """ Calculate the free-free Gaunt factor.
    @param Z : Ion charge

    @param electron_velocity : Velocity of the impinging electron
    @unit : Speed of light c

    @param photon_energy : The energy of the photon
    @unit : Hartree
    """

    if use_BvdH:
        return gauntFactorFF_BvdH(Z, electron_velocity, photon_energy)

    v = electron_velocity
    en = photon_energy

    #eta = Z/sqrt(E/Ry)  = Z/sqrt(0.5mv**2/Ry) = Z/sqrt(mv**2/Ha) = Z/c/sqrt(v**2/c**2) = Z*alpha/v
    eta1 = Z * alpha / v

    #E = hbar^2 k ^2 / 2m =(Ha) (k aB)^2/2
    #=> k = sqrt(2 E) = sqrt( 2* 0.5*mv^2 )  = sqrt(v^2/c^2)*c = v/alpha

    k1 = v / alpha
    # 0.5*mvf**2 = 0.5*mvi**2 - en
    #mvf^2 = mvi^2 + 2*en
    #vf^2 = vi^2 + 2en/m
    #vf = sqrt(vi^2 + 2en/m)
    #vf/c = sqrt(vi^2/c^2 + 2*en/c^2) = sqrt(vi^2 + 2*en*alpha^2)
    eta2 = Z * alpha / numpy.sqrt(v**2 + 2.*en * alpha**2)

    v_f = Z * alpha / eta2
    k2 = v_f / alpha

    prefactor = 2.*3**0.5/math.pi
    I0 = I(eta1, eta2, k1, k2, 0)
    I1 = I(eta1, eta2, k1, k2, 1)

    tmp = 1./eta1/ eta2 * ( (eta1**2 + eta2**2 + 2.*eta1**2 * eta2**2) * I0 - 2.*eta1 * eta2 * (1. + eta1**2)**0.5  * (1. + eta2**2)**0.5 * I1) * I0

    return prefactor * tmp.real

def tauMinus(n1, n2, l):
    if l == 0:
        return 0.0
    tmp = 2**(2.*l) / factorial(2.*l - 1)
    tmp *= ( gamma(n1+l+1) * gamma(n2+l) / gamma(n1 - l) / gamma(n2-l+1) )**0.5
    tmp *= (n1*n2)**(l+1)
    tmp /= (n1 + n2)**(n1+n2)
    tmp *= (n1 - n2)**(n1 - 2 - l)
    tmp *= (n2 - n1)**(n2 - l)

    tmp2 = hyp2f1(l+1-n1, l-n2, 2*l, -4.*n1*n2/(n1-n2)**2)
    tmp2 -= ((n1-n2)/(n1+n2))**2 * hyp2f1(l-1-n1, l-n2, 2*l, -4.*n1*n2/(n1-n2)**2)

    return tmp * tmp2

def tauPlus(n1, n2, l):
    return tauMinus(n2, n1, l+1)

def crossSectionBoundFree(Z, n, l, e_kin):
    prefactor = 16.*math.pi**2/3. * alpha**2 * 2. * e_kin
    n1 = n
    n2 = 1j * math.sqrt(2./e_kin)
    term1 = (l + 1.)/(2.*l + 1) * tauPlus(n1, n2, l)**2
    term2 = l/(2.*l + 1.) * tauMinus(n1, n2,l)**2

    return prefactor*( term1 + term2 )

def crossSectionBoundFreeClassical(Z, n, e_kin):
    v = math.sqrt(2.*e_kin)*alpha
    prefactor = 16./3./math.sqrt(3.) * alpha**2
    n1 = n
    n2 = 1j * math.sqrt(2./e_kin)
    rho = n2 / n1

    return prefactor / v / n * ( rho**2 / (1.+rho**2) )**2

def gauntFactorBoundFree(Z, n, l, e_kin):
    return crossSectionBoundFree(Z, n, l, e_kin) / crossSectionBoundFreeClassical(Z, n, e_kin)

def oscillatorStrengthFreeFree(eta1, eta2, occupation1):
    """ Free-Free oscillator strength, Menzel-Pekeris (1.22)
    """
    prefactor = 2**6/3./occupation1

    nominator = math.exp(-2.*math.pi * eta1)
    denominator = eta1**2 * eta2**2 * (1./eta1**2  - 1./eta2**2)**3 * (eta2 - eta1) * (1. - math.exp(-2.*math.pi * eta1)) * (1. - math.exp(-2.*math.pi * eta2))

    return prefactor * nominator / denominator * abs(DeltaHypergeometric(1j*eta1, 1j*eta2))


def oscillatorStrengthFreeBound(eta1, n2):
    """ Free-Bound oscillator strength, Menzel-Pekeris (1.21)
    """
    prefactor = 2**5/3./n2**2

    nominator = math.exp( -4. * eta1 * math.arctan(n2/eta1) )
    denominator = eta1**3 * n2**3 * (1./eta1**2  + 1./n2**2)**3.5 * (1. - math.exp(-2.*math.pi * eta1))

    return prefactor * nominator / denominator * abs(DeltaHypergeometric(1j*eta1, n2))


def oscillatorStrengthBoundBound(n1, n2, occupation2):
    """ Bound-bound oscillator strength, Menzel-Pekeris (1.15)
    """
    prefactor = 2**6/3./occupation2

    nominator = ((n1 - n2)/(n1 + n2))**(2*n1 + 2*n2)
    denominator = n1**2 * n2**2 * (1./n2**2 - 1./n1**2)**3 * (n1 - n2)

    return abs( prefactor * nominator / denominator * DeltaHypergeometric(n1, n2) )


def DeltaHypergeometric(x1, x2):
    """ Difference of squared hypergeometric functions evaluated at x1 and x2, respectively.

    @param x1 : First argument.
    @param x2 : Second argument.

    @return   : | F(-x1+1, -x2, 1, -4 x1 x2 /(x1-x2)^2) |^2 - | F(-x2+1, -x1, 1, -4 x1 x2 /(x1-x2)^2) |^2
    """

    # x variable:
    y = -4. * x1 * x2 / (x1 - x2)**2


    f1 = (hyp2f1(-x1+1, -x2, 1., y))**2
    f2 = (hyp2f1(-x2+1, -x1, 1., y))**2

    Delta = f1 - f2

    return Delta


#def gauntFactorFFAvg(






