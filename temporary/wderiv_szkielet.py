import wderiv
import numpy as np

# ver1 duża funckja
f = # funkcja 2D, np. fx=x**2*np.arrange(-n/2, n/2), fy =1/3*y**2*np.arrange(-n/2, n/2), f=np.meshgrid(fx,fy)
fd = wderiv.derivative(function=f, rank=2)
wderiv.plot(fd)

# ver2 klasa
fnp = ...
f = WDeriv(fnp)# funkcja 2D, stworzona przez lib wderiv [-L/2, L/2]
f.derivative(rank=2)
f.plot()

###########
# WDERIV
# To jest binding C do Pythona
#! Do przemyślenia, czy jest to jedna duża funkcja ver1 czy klasa ver2
def derivative(function, rank=1): #  function to funkcja wyjściowa z numpy #! ver1
    type_ = type(function) # rodzaj zmiennej, real albo complex
    dim_ = function.shape()
    if type_ == np.complex and dim_ = 2:
        wderiv_d2fdx2_2d_c(...) # wywołanie funkcji z C

class WDeriv: #! ver2
    self.function=...
    self.type_=type(self.function)
    self.dim_=len(self.function.shape)
    ...
    def derivative(self, ...):
        if self.type_ == np.complex and self.dim_ = 2:
            wderiv_d2fdx2_2d_c(...) # wywołanie funkcji z C
    ...
    def plot(self, ...):
        yrange=self.shape[1] #???
        plt.set_ylim