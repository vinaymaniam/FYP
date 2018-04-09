import numpy as np

from ufresh_helper_functions import heir_to_standard, heirarchical_search, reconstruct_from_map
from thirdparty import col_to_im, im_to_col, merge_patch

def ufresh(xtest, blocksize, heirarchy, index, mymap):
    cropwidth = [len(xtest),len(xtest[0])]

    xtestvec = im_to_col(xtest, blocksize, 1)
    dcx = np.mean(xtestvec, axis=0)
    xtestvec = xtestvec - np.tile(dcx, [xtestvec.shape[0], 1])

    ind = heirarchical_search(xtestvec, heirarchy)
    idx = heir_to_standard(ind, index)

    xrecmean = reconstruct_from_map(xtestvec, mymap, idx, dcx)

    xrecim = merge_patch(xrecmean, blocksize, cropwidth).T

    return xrecim

# a = np.random.rand(3,3)
# b = np.random.rand(3,3)
# print(a)
# print(b)
# print(a.mean(axis=0).mean(axis=0))
# mse = abs(a-b).mean(axis=0).mean(axis=0)
# psnr = 10*np.log10(1/mse)
# print(psnr)

# a = np.identity(3)
# print(a)
# print(a*3)