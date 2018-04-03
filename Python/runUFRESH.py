import numpy as np
import time
import scipy as sp

import pywt

from utils import myglob, load_images
from backprojection_2X import backprojection_2X
from ufresh import ufresh

directory_x = 'Testing_Images/FRESH_upscaled/Set5'
pattern = '*.bmp'
directory_y = 'Testing_Images/GT/Set5'

XpathCell = myglob(directory_x, pattern)
Xcell = load_images(XpathCell)
YpathCell = myglob(directory_y, pattern)
Ycell = load_images(YpathCell)

blocksize = [5, 5]
stepsize = [1, 1]

Psnr = np.ndarray([1, len(Xcell)])

for imgIdx in range(0,len(Xcell),1):
    t1 = time.time()
    print('--------------------------------------------------------')
    print('Processing image ', str(imgIdx), ' of total ', str(len(Xcell)))
    Xtest = Xcell[imgIdx]
    Ytest = Ycell[imgIdx]

    for stage in [1,2]:
        heir = sp.io.loadmat('pyHeirarchy4096.mat')
        heirarchy = heir['heirarchy']
        index = heir['index']
        mymap = sp.io.loadmat('pyMap4096cell96.mat')
        Map = mymap['Map']

        Xrec = np.zeros([len(Xtest), len(Xtest[0], 4)])
        for rot in range(0,4):
            Xtestrot = sp.ndimage.interpolation.rotate(Xtest, 90*rot)
            X = ufresh(Xtestrot, blocksize, stepsize, heirarchy, index, Map)
            X = sp.ndimage.interpolation.rotate(X, 360 - 90*rot)
            X = backprojection_2X(X, Ytest)

        ensembleMean = np.mean(Xrec,axis=2)
        Xtest = backprojection_2X(ensembleMean, Ytest)

    mse = abs(Xtest-Ytest).mean(axis=0).mean(axis=0)
    Psnr[imgIdx] = 10*np.log10(1/mse)
    print(Psnr[imgIdx])
    t2 = time.time()
    print('Elapsed time is ', str(t2-t1))

print('Average PSNR across all runs = ', str(Psnr.mean(axis=0)))






