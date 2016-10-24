import matplotlib
import matplotlib.pyplot as plt
import numpy

data100 = numpy.loadtxt('data_100percent.dat')
l100 = data100[:,0]
data100 = data100[:,1]
fit100 = numpy.loadtxt('fit_100percent.dat')[:,1]

data100_1 = numpy.loadtxt('data_100percent_T1.dat')
l100_1 = data100_1[:,0]
data100_1 = data100_1[:,1]
fit100_1 = numpy.loadtxt('fit_100percent_T1.dat')[:,1]

data30 = numpy.loadtxt('data_30percent.dat')
l30 = data30[:,0]
data30 = data30[:,1]
fit30 = numpy.loadtxt('fit_30percent.dat')[:,1]

data10 = numpy.loadtxt('data_10percent.dat')
l10 = data10[:,0]
data10 = data10[:,1]
fit10 = numpy.loadtxt('fit_10percent.dat')[:,1]

# Offsets.
data30 += 0.4
fit30 += 0.4
data100 += 0.6
fit100 += 0.6
fit100_1 += 0.6

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)

linewidth = 2.0
plt.plot(l10, data10, 'k.')
plt.plot(l10, fit10, 'k', linewidth=linewidth)
plt.plot(l30, data30, 'g.')
plt.plot(l30, fit30, 'k', linewidth=linewidth)
plt.plot(l100, data100, 'r.')
plt.plot(l100, fit100, 'k', linewidth=linewidth)
plt.plot(l100, fit100_1, 'k--', linewidth=linewidth)

plt.xlabel(r'photon wavelength (nm)')
plt.ylabel(r'photon flux (a.u.)')

plt.ylim([-0.2, 1.0])

plt.text( x=14.0, y=0.0, s='power = 10 %')
plt.text( x=13.0, y=0.4, s='power = 30 %', color='g')
plt.text( x=12.0, y=0.8, s='power = 100 %', color='r')

plt.show()


