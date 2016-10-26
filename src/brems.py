#! /usr/bin/env python2.7

import numpy, math
import argparse
from scipy import constants
from scipy.constants import physical_constants as PC
from matplotlib import pyplot
from GauntCalculator import freeFreeThermalSpectrum

Ry_in_J = constants.value('Rydberg constant times hc in J')
Hz_in_eV = constants.h/constants.e
aB_in_cm = constants.value('Bohr radius') * 1.e2

def main():
    """ The main routine."""

    parser = argparse.ArgumentParser(description='Calculate thermal bremsstrahlung emission.')

    # Ion density.
    parser.add_argument('ne', action='store', type=float, help="The electron density in cm^-3")
    # Average charge
    parser.add_argument('Z', action='store', type=float, help="The average ionization")
    # Temperature
    parser.add_argument('T', action='store', type=float, help="The plasma temperature in eV")
    # Spectrum or total
    parser.add_argument('-s','--spectrum',  action='store', default='energy', help="Which spectrum to compute",  dest='calculate_spectrum')
    #parser.add_argument('-t','--total',  action='store_true', default=True, help="Whether to compute the total emitted power", dest='calculate_total')
    parser.add_argument('-p','--plot',  action='store_true', default=False, help="Whether to plot the spectrum", dest='plot')

    # Parse command line arguments.
    args = parser.parse_args()

    print ""
    print "********************************************************************************"
    print "* Starting bremsstrahlung calculation for"
    print "* ne = %4.3e/cm^3 " % (args.ne)
    print "* Z = %4.3f " % (args.Z)
    print "* T = %4.3f eV " % (args.T)
    print ""

    # Setup energy range. We want the spectrum to go from ~10^-2 to 10^2 T
    lg_T = math.log(args.T, 10)
    lg_T_lo = math.floor(lg_T) - 2
    lg_T_hi = math.floor(lg_T) + 2

    Ha_to_eV = 27.2
    lg_energies_eV = numpy.linspace(lg_T_lo, lg_T_hi, 41)
    energies_eV = 10.**lg_energies_eV

    print "* Setting energy range to %4.3f to %4.3f eV using logarithmic binning." % (energies_eV[0], energies_eV[-1])

    energies_Ha = energies_eV/Ha_to_eV

    ne_aB3 = args.ne * aB_in_cm**3
    ni_aB3 = ne_aB3 / args.Z

    print "* Calculating free-free thermal emisssion ..."
    spectrum = freeFreeThermalSpectrum(args.Z, args.T/Ha_to_eV, energies_Ha) * ne_aB3 * ni_aB3 # Ry/s/Hz/cm^3
    # Convert to W/eV/cm^3
    spectrum = spectrum * Ry_in_J / Hz_in_eV / aB_in_cm**3
    # Convert to eV/s/Hz
    #spectrum = spectrum * 13.6 / Hz_in_eV
    print "* done."
    ordinate = energies_eV
    abcissa = spectrum

    ordinate_title = "Photon energy [eV]"
    abcissa_title = "Power spectrum dP/dE dV [W/eV/cm^3]"

    if args.calculate_spectrum.lower() == "wavelength":
        # Convert to wavelength spectrum.
        print "* Converting to wavelength spectrum."
        eV_to_nm = 1239.8
        wavelengths_nm = eV_to_nm / energies_eV
        spectrum = eV_to_nm * spectrum / wavelengths_nm**2

        ordinate = wavelengths_nm[::-1]
        abcissa = spectrum[::-1]

        ordinate_title = "Wavelength [nm]"
        abcissa_title = "Power spectrum dP/dlambda dV [W/nm/cm^3]"

    outdata = zip(ordinate, abcissa)

    ##print outdata
    header="Bremsstrahlung spectrum\n\
ne = %4.3e/cm^3\n\
Z = %4.3f\n\
T = %4.3f eV\n\
\n\
%s %s\n" % (args.ne, args.Z, args.T, ordinate_title, abcissa_title) % ()

    filename = 'brems_ne%4.3e_Z%4.3f_T%4.3f.txt' % (args.ne, args.Z, args.T)
    print "* Saving data to %s." % (filename)
    numpy.savetxt(filename, outdata, header=header, fmt='%.6e' )

    if args.plot:
        print "* Setting up plot."
        pyplot.loglog( ordinate, spectrum )
        pyplot.xlabel( ordinate_title )
        pyplot.ylabel( abcissa_title )
        pyplot.title("Bremsstrahlung spectrum \n\
ne = %4.3e/cm^3, Z = %4.3f, T = %4.3f eV" % (args.ne, args.Z, args.T) )

        pyplot.show()

    print "* All done."
    print "********************************************************************************"

if __name__ == "__main__":
    main()
