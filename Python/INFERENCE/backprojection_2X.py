import numpy as np
import scipy as sp
import pywt

# DOESNT WORK FOR RECTANGULAR IMAGES!!!!!!
def backprojection_2X(Ir, Orig):
    myfilter = 'db2'#'bior4.4'
    upscale_lvl = 1

    mode = 'smooth'
    wd1 = pywt.wavedec2(Orig,myfilter,level=upscale_lvl, mode=mode)
    # corg = np.append(wd1[0].ravel(),[wd1[1][0].ravel(), wd1[1][1].ravel(), wd1[1][2].ravel()])
    ilowc = wd1[0]/2
    rangeImg = [np.amin(ilowc), np.amax(ilowc)]

    Ir = range0toN(Ir, rangeImg)

    wd2 = pywt.wavedec2(Ir, myfilter, level=upscale_lvl, mode=mode)

    crec = wd2
    crec[0] = 2*ilowc
    if(wd2[0].shape != wd2[1][0].shape):
        lst = list()
        for i in range(0,len(wd2[1])):
            # lst.append(wd2[1][i].T)
            lst.append(np.reshape(wd2[1][i].ravel(), wd2[0].shape))
        crec[1] = lst


    irec = pywt.waverec2(crec, myfilter, mode=mode)
    ibp = range0toN(irec, rangeImg)
    return ibp

def range0toN(A, myrange):
    a = myrange[0]
    b = myrange[1]
    # A[:] = [b if ele > b else a if ele < a else ele for ele in A]
    A[A > b] = b
    A[A < a] = a
    return A