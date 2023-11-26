import numpy as np
import cupy as cp
import time
# Matrix multiplication with NumPy
start_time = time.time()
np_array = np.random.rand(100,100,100,6,6)
np_result = np.matmul(np_array, np_array)
np_time = time.time() - start_time
# Matrix multiplication with CuPy
start_time = time.time()
cp_array = cp.random.rand(100,100,100,6,6)
cp_result = cp.matmul(cp_array, cp_array)
cp_time = time.time() - start_time
print(f"NumPy Time: {np_time}")
print(f"CuPy Time: {cp_time}")