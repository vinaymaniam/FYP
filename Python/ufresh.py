import numpy as np

from ufresh_helper_functions import heir_to_standard, heirarchical_search, reconstruct_from_map

def ufresh():
    return 0

# a = np.random.rand(3,3)
# b = np.random.rand(3,3)
# print(a)
# print(b)
# print(a.mean(axis=0).mean(axis=0))
# mse = abs(a-b).mean(axis=0).mean(axis=0)
# psnr = 10*np.log10(1/mse)
# print(psnr)

a = np.identity(3)
print(a)
print(a*3)