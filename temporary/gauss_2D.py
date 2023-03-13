'''
Developed within project Spectral Derivers

Koło naukowe fizyki komputerowej
Authors: Viktoriia Vlasenko
'''


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def Gauss2D(x=0, y=0, mx=0, my=0, sx=1, sy=1):
	return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))

def GaussDerivative_X(x=0, y=0, mx=0, my=0, sx=1, sy=1):
	return 1. / (2. * np.pi * sx**3 * sy) * (-x + mx)  * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))

def GaussDerivative_Y(x=0, y=0, mx=0, my=0, sx=1, sy=1):
	return 1. / (2. * np.pi * sx * sy**3) * (-y + my)  * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))	


if __name__ == '__main__':
	sigma  = 1
	mean   = 0
	n = 500
	L = 30
	dx = L/n
	x = np.arange(-L/2, L/2, dx, dtype = 'complex_')
	y = np.arange(-L/2, L/2, dx, dtype = 'complex_')
	x, y = np.meshgrid(x, y)
	z = Gauss2D(x, y)

	df_x = GaussDerivative_X(x, y)
	df_y = GaussDerivative_Y(x, y)

	fhat = np.fft.fft2(z)
	kappa_x = (2*np.pi/L)*np.arange(-n/2, n/2)
	kappa_y = (2*np.pi/L)*np.arange(-n/2, n/2)
	kappa_x, kappa_y = np.meshgrid(kappa_x, kappa_y)
	kappa_x = np.fft.fftshift(kappa_x)
	dfhat = kappa_x*fhat*(1j)
	dfFFT = np.real(np.fft.ifft2(dfhat))
	fig = plt.figure()
	ax1 = plt.subplot(2, 2, 1)
	ax1.contourf(x, y, z, cmap=cm.viridis)
	ax1.set_xticks([])
	ax1.set_yticks([])
	ax1.set_title('2D Gauss function')
	ax1.set_xlabel(r'$x_1$')
	ax1.set_ylabel(r'$y_1$')

	ax2 = plt.subplot(2, 2, 2)
	ax2.contourf(x, y, df_x, cmap=cm.viridis)
	ax2.grid(False)
	ax2.set_xticks([])
	ax2.set_yticks([])
	ax2.set_title(r'$\frac{\delta F}{\delta x}$')
	ax2.set_xlabel(r'$x_1$')
	ax2.set_ylabel(r'$y_1$')

	ax3 = plt.subplot(2, 2, 3)
	ax3.contourf(x, y, dfFFT, cmap=cm.viridis)
	ax3.grid(False)
	ax3.set_xticks([])
	ax3.set_yticks([])
	ax3.set_title('FFT derivative')
	ax3.set_xlabel(r'$x_1$')
	ax3.set_ylabel(r'$y_1$')
	# ax4 = plt.subplot(2, 2, 4)
# 	# ax[0].plot_surface(x, y, z, rstride=3, cstride=3, linewidth=1, antialiased=True,
#  #                cmap=cm.viridis)
	
	plt.savefig(f"gauss2D.jpg", dpi=150)