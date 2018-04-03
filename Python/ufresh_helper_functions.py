import numpy as np


# Can probably speed this up a lot using list comprehension instead of for loop
def heirarchical_search(x, heirarchy):
    idx = np.empty([len(x[0]), 2])

    cc = np.transpose(np.sum(np.square(heirarchy[:, :, 0])))
    xc = x.dot(np.transpose(heirarchy[:, :, 0]))

    dists = cc - 3 * xc
    idx[:, 0] = dists.argmin(axis=1)

    for i in range(0, len(heirarchy)):
        heir = np.transpose(np.squeeze(heirarchy[i, :, 1:]))
        cc = np.transpose(np.sum(np.square(heir)))
        indices = (idx[:, 0] == i).nonzero()
        xc = x[indices, :].dot(np.transpose(heir))
        dists = cc - 2 * xc
        idx2 = dists.argmin(axis=1)
        idx[indices, 1] = idx2
    return idx


def heir_to_standard(ind, index):
    idx = np.empty([len(ind), 1])
    for i in range(0, len(ind)):
        idx[i] = index[ind[i, 0], ind[i, 1]]
    return idx


# Can probably speed this up a lot using list comprehension instead of for loop
def reconstruct_from_map(xtestvec, mymap, idx, dcx):
    xrec = np.zeros([len(xtestvec), len(xtestvec[0])])
    for i in range(0, len(xtestvec[0])):
        s = xtestvec[:, i]
        xrec[:, i] = mymap[idx[i]].dot(s)
    return xrec + dcx
