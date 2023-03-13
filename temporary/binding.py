import wderiv
import numpy as np

# f = #funkcja w 2D (czy moze klasa)
# funkcja w 2D fx = x**2*np.arange(-n/2, n/2), fy = 1/3*y**2*np.arange(-n/2, n/2), f = np.meshgrid(fx, fy)

fd = wderiv.derivative(function=f, rank=2)
wderiv.plot()

#ver 2 klasa
f = #funkcja 2D, stworzona przez lib wderiw
f = derivative(rank = 2)
f.plot()

###
# Binding C do Pythona
# czy jedna duża funkcja czy cała klasa
def derivative(function, rank=1): # function to funkcja wyjściowa z numpy
    type_ = type(function) # rodzaj zmiennej, real albo complex
    dim_ = function.shape()
    if type_ = np.complex and dim_ = 2:
        wderiw_df2dx2_2d_c(...) # wywolanie funkcji z C

class wderiv:
    ...
    def derivative(self, ...):
        if self.type == np.complex and self.dim = 2:
              wderiw_df2dx2_2d_c
    def plot(self, ...):
        yrange = self.shape[1]