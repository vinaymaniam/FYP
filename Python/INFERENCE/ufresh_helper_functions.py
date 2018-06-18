import numpy as np
from numpy import *

# Can probably speed this up a lot using list comprehension instead of for loop
def heirarchical_search(x, heirarchy):
    idx = np.empty([len(x[0]), 2], dtype=int)
    x = x.T
    x = x.astype(float32)

    cc = np.transpose(np.sum(np.square(heirarchy[:, :, 0])))
    xc = x.dot(np.transpose(heirarchy[:, :, 0]))

    dists = cc - 2 * xc
    idx[:, 0] = dists.argmin(axis=1)

    for i in range(0, len(heirarchy)):
        heir = (np.squeeze(heirarchy[i, :, 1:])).T
        cc = (np.sum(np.square(heir), axis=1)).T
        indices = (idx[:, 0] == i).nonzero()
        # xc = np.squeeze(x[indices, :].dot(heir.T))
        xc = x[indices, :].dot(heir.T)
        dists = np.reshape(cc - 2 * xc,[-1,len(cc)])
        idx2 = dists.argmin(axis=1)
        idx[indices, 1] = idx2
    return idx


def heir_to_standard(ind, index):
    idx = np.empty([len(ind), 1], dtype=int)
    for i in range(0, len(ind)):
        idx[i] = index[ind[i, 0], ind[i, 1]]
    return idx


# Can probably speed this up a lot using list comprehension instead of for loop
def reconstruct_from_map(xtestvec, mymap, idx, dcx):
    xtestvec = xtestvec.astype(float32)
    xrec = np.zeros([len(xtestvec), len(xtestvec[0])],dtype=float32)
    for i in range(0, len(xtestvec[0])):
        s = xtestvec[:, i]
        xrec[:, i] = mymap[int(idx[i])].dot(s)
    return xrec + np.tile(dcx, [xrec.shape[0], 1])
