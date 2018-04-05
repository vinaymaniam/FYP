import numpy as np


def im_to_col(A, blocksize, stepsize=1):
    # Parameters
    m, n = A.shape
    s0, s1 = A.strides
    nrows = m - blocksize[0] + 1
    ncols = n - blocksize[1] + 1
    shp = blocksize[0], blocksize[1], nrows, ncols
    strd = s0, s1, s0, s1

    out_view = np.lib.stride_tricks.as_strided(A, shape=shp, strides=strd)
    return out_view.reshape(blocksize[0] * blocksize[1], -1)[:, ::stepsize]


def col_to_im(B, block_size, image_size):
    m, n = block_size
    mm, nn = image_size
    return B.reshape([nn - n + 1, mm - m + 1]).T

def merge_patch(p, blocksize, cropwidth):
    y = blocksize[0]
    x = blocksize[1]
    Y = cropwidth[0]
    X = cropwidth[1]
    img = np.zeros([Y,X])
    coeff = np.zeros([Y,X])
    p_idx = 0

    for xx in range(0,x):
        for yy in range(0,y):
            pp = col_to_im(p[p_idx,:],[y,x],[Y,X])
            img[yy:yy+Y-y+1,xx:xx+X-x+1] = img[yy:yy+Y-y+1,xx:xx+X-x+1] + pp
            coeff[yy:yy+Y-y+1,xx:xx+X-x+1] = coeff[yy:yy+Y-y+1,xx:xx+X-x+1] + 1
            p_idx = p_idx + 1
    img = np.divide(img,coeff)
    return img


def count_cover(sz, blocksize, stepsize):
    cnt = np.ones(sz)
    for k in range(0, len(sz)):
        ids = np.arange(sz[k]).transpose()
        s = sz
        s = np.delete(sz,k)
        ids = np.reshape(ids[:, np.ones([1, np.prod(s)])], [len(ids), s])

        cnt = cnt * np.maximum(np.min(np.floor((ids-1)/stepsize[k]),
                            np.floor((sz[k]-blocksize[k]))) -
                            np.maximum(np.ceil((ids-blocksize[k])/stepsize[k]),0)+1, 0)
