'''
Developed within project Spectral Derivers

Koło naukowe fizyki komputerowej
Authors: Viktoriia Vlasenko, Wojciech Grunwald
'''


import matplotlib.pyplot as plt
import numpy as np
import scipy as sci


def Gauss(x, mean, sigma):
	return np.exp(-(x - mean)**2 / (2 * sigma**2))/sigma/((sci.pi*2)**(-1/2))

def GaussDerivative(x, mean = 1, sigma = 1):
	return np.exp(-(x - mean)**2 / (2 * sigma**2))/sigma/((sci.pi*2)**(-1/2))*(mean-x)/sigma**2

def d_Gauss(x, mean, sigma):
	# numerical derivative
    return (Gauss(x+dx, mean, sigma) - Gauss(x, mean, sigma))/dx

def Fft(f):
    return np.fft.fftshift(np.abs(np.fft.fft(f)))

def Ifft(fftf):
    return np.abs(np.fft.fftshift(np.fft.ifft(fftf)).real)


if __name__ == '__main__':
	sigma  = 3/2
	mean   = 0
	X = np.arange(-10, 10, 0.01)
	dx = np.abs(X[1] - X[0])
	# freq_x = np.fft.fftshift(np.fft.fftfreq(np.shape(X)[0])) / dx * 2 * sci.pi

	Y1 = [Gauss(x, mean, sigma) for x in X]
	Y2 = [GaussDerivative(x, mean, sigma) for x in X]
	Y3 = [d_Gauss(x, mean, sigma) for x in X]
	Y4 = Fft(Y3)
	Y4 = Y4*(1j)
	Y4 = Ifft(Y3)

	fig, ax = plt.subplots(2, figsize=(8, 6))
	ax[0].scatter(X, Y1, marker = "o", s=2, color='orangered', label = 'Gauss function')
	ax[0].scatter(X, Y2, marker = "o", s=2, color='green', label = 'Analytic derivative')
	ax[1].scatter(X, Y3, marker = "o", s=2, color='red', label = 'Numerical derivative')
	ax[1].scatter(X, Y4, marker = "o", s=2, color='blue', label = 'FTT derivative')

	ax[0].set_xlabel("X")
	ax[0].set_ylabel("Y", rotation=0, labelpad=20)
	ax[0].minorticks_on()
	ax[0].grid(which='major')
	ax[0].grid(which='minor', linestyle=':')
	ax[0].legend()

	ax[1].set_xlabel("X")
	ax[1].set_ylabel("Y", rotation=0, labelpad=20)
	ax[1].minorticks_on()
	ax[1].grid(which='major')
	ax[1].grid(which='minor', linestyle=':')
	ax[1].legend()
	plt.show()
