import numpy as np
import time
import scipy as sp
from scipy import io, ndimage
from numpy import *


import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pywt

from utils import myglob, load_images
from backprojection_2X import backprojection_2X
from ufresh import ufresh


def runUFRESH():
    directory_x = 'Testing_Images/FRESH_upscaled/Set5'
    pattern = '.bmp'
    directory_y = 'Testing_Images/GT/Set5'

    XpathCell = myglob(directory_x, pattern)
    Xcell = load_images(XpathCell)
    YpathCell = myglob(directory_y, pattern)
    Ycell = load_images(YpathCell)

    blocksize = [5, 5]

    Psnr = np.zeros([len(Xcell)])
    prepsnr = np.zeros([len(Xcell)])
    for imgIdx in range(0,len(Xcell)):
        t1 = time.time()
        print('--------------------------------------------------------')
        print('Processing image ', str(imgIdx+1), ' of total ', str(len(Xcell)))
        Xtest = Xcell[imgIdx]
        Ytest = Ycell[imgIdx]
        mse = np.square((Xtest - Ytest).ravel()).mean(axis=0)
        prepsnr[imgIdx] = 10 * np.log10(1 / mse)
        print('PSNR before processing = ', prepsnr[imgIdx])
        for stage in [1,2]:
            heir = sp.io.loadmat('pyHeirarchy4096.mat')
            heirarchy = heir['heirarchy'].astype(float32)
            index = heir['index'].astype(float32)
            mymap = sp.io.loadmat('pyMap4096cell96.mat')

            C = np.squeeze(mymap['Map'])
            # Map = np.empty((C.shape[0], C[0].shape[0], C[0].shape[1]))
            # for i in range(Map.shape[0]):
            Map = [np.ndarray([C[0].shape[0], C[0].shape[1]])]*C.shape[0]
            for i in range(len(Map)):
                Map[i] = C[i].astype(float)

            Xrec = np.zeros([len(Xtest), len(Xtest[0]), 4])
            for rot in range(0,4):
                print(rot)
                Xtestrot = np.rot90(Xtest, rot)
                # Xtestrot = sp.ndimage.rotate(Xtest, 90*rot)
                X = ufresh(Xtestrot, blocksize, heirarchy, index, Map)
                X = backprojection_2X(X, np.rot90(Ytest, rot))
                X = np.rot90(X, 4-rot)
                # X = sp.ndimage.rotate(X, 360-90*rot)
                Xrec[:,:,rot] = X

            Xtest = np.mean(Xrec,axis=2)
            Xtest = backprojection_2X(Xtest, Ytest)

        mse = np.square((Xtest-Ytest).ravel()).mean(axis=0)
        Psnr[imgIdx] = 10*np.log10(1/mse)
        print('PSNR after processing  = ',Psnr[imgIdx])
        t2 = time.time()
        print('Elapsed time is ', str(t2-t1))

        # plt.figure()
        # plt.imshow(Xtest, cmap='gray')
        # plt.title('Xtest')
        # plt.figure()
        # plt.imshow(Ytest, cmap='gray')
        # plt.title('Ytest')
        # plt.show()
    print('========================================================')
    print('Average PSNR across all runs = ', str(Psnr.mean(axis=0)))
    print('Average improvement in PSNR  = ', str((Psnr-prepsnr).mean(axis=0)))






