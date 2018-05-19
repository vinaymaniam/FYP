import numpy as np
from numpy.linalg import inv
import scipy
import scipy.io
import random
from tqdm import tqdm



def mapping_calculation(dxloc,dyloc,i,clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    X = xdict['X']
    Y = ydict['Y']
    loadfilename = 'data_files/pyCenter'+str(i)+'.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    # clusterszA = 96-1
    lam = 0.01
    Map = np.ndarray([cn,25,25])
    Res = np.ndarray([cn,25])
    for t in tqdm(range(0,cn,1)):
        c1 = Center[...,t].reshape(1,-1)
        # This cdist part takes a very long time
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        LR = np.squeeze(X[..., idx[0:clusterszA]], axis=2)
        HR = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
        M = HR.dot(LR.transpose().dot(inv(LR.dot(LR.transpose())
                   + lam * np.identity(len(Center)))))
        Map[t,:,:] = M
        Res[t,:] = HR.mean(axis=1) - np.dot(M,LR.mean(axis=1))
    savfilename = 'data_files/pyMap'+str(i)+'mat'+str(clusterszA+1)+'.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return

def mapping_calculationRANSAC(dxloc,dyloc,i,clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    # Mean removed versions for mapping calculation
    X = xdict['X']
    Y = ydict['Y']
    # Non mean removed versions for PSNR calculation in RANSAC
    X1 = xdict['X1']
    Y1 = ydict['Y1']

    loadfilename = 'data_files/pyCenter'+str(i)+'.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    # clusterszA = 96-1
    lam = 0.01
    Map = np.ndarray([cn,25,25])
    Res = np.ndarray([cn,25])
    for t in tqdm(range(0,cn,1)):
        c1 = Center[...,t].reshape(1,-1)
        # This part takes a most of the time so the RANSAC shouldn't add too much to that
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        # DO RANSAC HERE----------------------------------------------------------------------------------------
        indices = idx[0:clusterszA]
        numTrials = 200 # This bit takes 0.6% of the total time so 200 trials wont severely affect training time
        nRansac = 10     # Size of subset to take from indices for each iteration of ransac
        thRansac = 1.2 # Threshold for detecting outliers(PSNR_best/PSNR_current > thRansac then outlier)
        minoutliers = np.inf
        Mbest = np.zeros([len(Center), len(Center)])
        for trial in range(0,numTrials):
            subsample = random.sample(range(0,len(indices)),nRansac)
            ind = indices[subsample]
            LR = np.squeeze(X[..., ind], axis=2)
            HR = np.squeeze(Y[..., ind], axis=2)
            M = HR.dot(LR.transpose().dot(inv(LR.dot(LR.transpose())
                       + lam * np.identity(len(Center)))))
            LR1 = np.squeeze(X1[..., indices], axis=2)
            HR1 = np.squeeze(Y1[..., indices], axis=2)
            SR1 = np.dot(M,LR1)

            mse = ((HR1 - SR1) ** 2).mean(axis=0)
            psnr = 10*np.log10(mse**-1)
            maxpsnr = psnr.max()
            noutliers = sum((maxpsnr/i) > thRansac for i in psnr)
            if(noutliers < minoutliers):
                minoutliers = noutliers
                Mbest = M
        print("Min Outliers = {}".format(minoutliers))#, end="\r")
        # END RANSAC HERE---------------------------------------------------------------------------------------
        Map[t,:,:] = Mbest
        LR = np.squeeze(X[..., indices], axis=2)
        HR = np.squeeze(Y[..., indices], axis=2)
        Res[t,:] = HR.mean(axis=1) - np.dot(Mbest,LR.mean(axis=1))
    savfilename = 'data_files/pyMap'+str(i)+'mat'+str(clusterszA+1)+'.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return
