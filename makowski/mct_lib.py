from numpy.ctypeslib import ndpointer
from ctypes import *

libCalc = CDLL("./libmct.so")

# call C function to check connection
libCalc.connect()

# Caculate Vectors
calculate_vectors = libCalc.calculate_vectors
calculate_vectors.restype = None
calculate_vectors.argtypes = [
    ndpointer(c_double),
    ndpointer(c_double),
    ndpointer(c_double),
    ndpointer(c_double),
    ndpointer(c_double),
    c_long,
    ndpointer(c_double),
    ndpointer(c_double),
    ndpointer(c_double),
    ndpointer(c_double),
    ndpointer(c_double),
]
print("# Loaded interpreter of `calculate_vectors()`")

# Caculate Means
calculate_mean = libCalc.calculate_mean
calculate_mean.restype = c_double
calculate_mean.argtypes = [ndpointer(c_double), c_long]
print("# Loaded interpreter of `calculate_mean()`")

# Caculate Covariance
calculate_covariance = libCalc.calculate_covariance
calculate_covariance.restype = c_double
calculate_covariance.argtypes = [
    ndpointer(c_double),
    ndpointer(c_double),
    c_long,
    c_double,
    c_double,
]
print("# Loaded interpreter of `calculate_covariance()`")
