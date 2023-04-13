import timeit
import numpy as np

import mct_lib

NELEMENTS = 134217728

var = {}
# fileName="/home2/archive/MCT-2021/lab1/var1.dat"
for i in range(1, 3, 1):
    fileName = f"./var{i}.dat"
    var[i] = np.fromfile(fileName, dtype=np.double)
    print(f"# Vector {i} is loaded")

print(
    mct_lib.calculate_covariance(
        var[1],
        var[2],
        NELEMENTS,
        mct_lib.calculate_mean(var[1], NELEMENTS),
        mct_lib.calculate_mean(var[2], NELEMENTS),
    )
)

