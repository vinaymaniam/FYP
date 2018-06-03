import numpy as np
from numpy.linalg import inv
import scipy
import scipy.io
import random

from tqdm import tqdm
import matplotlib.pyplot as plt


def mapping_calculation_6by6(dxloc, dyloc, i, clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    X = xdict['X']
    Y = ydict['Y']
    # Non mean removed versions for PSNR calculation
    X1 = xdict['X1']
    Y1 = ydict['Y1']
    loadfilename = 'data_files/pyCenter' + str(i) + '.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    Map = np.ndarray([cn, 36, 37])
    Res = np.ndarray([cn, 36])
    lam = 0.0003
    for t in tqdm(range(0, cn, 1)):
        c1 = Center[..., t].reshape(1, -1)
        # This cdist part takes a very long time
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        LR = np.squeeze(X[..., idx[0:clusterszA]], axis=2)
        HR = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
        # Add in bias term so regression model learns bias
        LR = np.append(LR, np.ones([1,LR.shape[1]]), axis=0)
        LLT = LR.dot(LR.transpose())
        LLT = LLT + lam * np.identity(len(LR))
        eigvals = np.linalg.eig(LLT)
        if (any(eigvals[0] < lam)):
            LLT = LLT + lam*np.identity(len(LR))
        # else:
        #     print('All good')
        M = HR.dot(LR.transpose().dot(inv(LLT)))
        Map[t, :, :] = M
    savfilename = 'data_files/pyMap' + str(i) + 'mat' + str(clusterszA + 1) + '.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return

def mapping_calculation_3by3(dxloc, dyloc, i, clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    X = xdict['X']
    Y = ydict['Y']
    # Non mean removed versions for PSNR calculation
    X1 = xdict['X1']
    Y1 = ydict['Y1']
    loadfilename = 'data_files/pyCenter' + str(i) + '.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    Map = np.ndarray([cn, 9, 10])
    Res = np.ndarray([cn, 9])
    lam = 0.0003
    for t in tqdm(range(0, cn, 1)):
        c1 = Center[..., t].reshape(1, -1)
        # This cdist part takes a very long time
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        LR = np.squeeze(X[..., idx[0:clusterszA]], axis=2)
        HR = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
        # Add in bias term so regression model learns bias
        LR = np.append(LR, np.ones([1,LR.shape[1]]), axis=0)
        LLT = LR.dot(LR.transpose())
        LLT = LLT + lam * np.identity(len(LR))
        eigvals = np.linalg.eig(LLT)
        if (any(eigvals[0] < lam)):
            LLT = LLT + lam*np.identity(len(LR))
        # else:
        #     print('All good')
        M = HR.dot(LR.transpose().dot(inv(LLT)))
        Map[t, :, :] = M
    savfilename = 'data_files/pyMap' + str(i) + 'mat' + str(clusterszA + 1) + '.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return

def mapping_calculation(dxloc, dyloc, i, clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    X = xdict['X']
    Y = ydict['Y']
    # Non mean removed versions for PSNR calculation
    X1 = xdict['X1']
    Y1 = ydict['Y1']
    loadfilename = 'data_files/pyCenter' + str(i) + '.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    Map = np.ndarray([cn, 25, 26])
    Res = np.ndarray([cn, 25])
    lam = 0.0003
    for t in tqdm(range(0, cn, 1)):
        c1 = Center[..., t].reshape(1, -1)
        # This cdist part takes a very long time
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        LR = np.squeeze(X[..., idx[0:clusterszA]], axis=2)
        HR = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
        # Add in bias term so regression model learns bias
        LR = np.append(LR, np.ones([1,LR.shape[1]]), axis=0)
        LLT = LR.dot(LR.transpose())
        LLT = LLT + lam * np.identity(len(LR))
        eigvals = np.linalg.eig(LLT)
        if (any(eigvals[0] < lam)):
            LLT = LLT + lam*np.identity(len(LR))
        # else:
        #     print('All good')
        M = HR.dot(LR.transpose().dot(inv(LLT)))
        Map[t, :, :] = M
    savfilename = 'data_files/pyMap' + str(i) + 'mat' + str(clusterszA + 1) + '.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return

# WRONG
def mapping_calculation_quad(dxloc, dyloc, i, clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    X = xdict['X']
    Y = ydict['Y']
    loadfilename = 'data_files/pyCenter' + str(i) + '.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    Map = np.ndarray([cn, 25, 51])
    Res = np.ndarray([cn, 25])
    lam = 0.0003
    Xq = np.append(X,np.square(X),axis=0)
    for t in tqdm(range(0, cn, 1)):
        c1 = Center[..., t].reshape(1, -1)
        # This cdist part takes a very long time
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        LR = np.squeeze(Xq[..., idx[0:clusterszA]], axis=2)
        HR = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
        # Add in bias term so regression model learns bias
        LR = np.append(LR, np.ones([1,LR.shape[1]]), axis=0)
        LLT = LR.dot(LR.transpose())
        eigvals = np.linalg.eig(LLT)
        if (any(eigvals[0] < lam)):
            LLT = LLT + lam*np.identity(len(LR))
        else:
            print('All good')
        M = HR.dot(LR.transpose().dot(inv(LLT)))
        Map[t, :, :] = M
    savfilename = 'data_files/pyMap' + str(i) + 'mat' + str(clusterszA + 1) + '.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return

def mapping_calculationRANSAC(dxloc, dyloc, i, clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    # Mean removed versions for mapping calculation
    X = xdict['X']
    Y = ydict['Y']
    # Non mean removed versions for PSNR calculation in RANSAC
    X1 = xdict['X1']
    Y1 = ydict['Y1']

    loadfilename = 'data_files/pyCenter' + str(i) + '.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    # clusterszA = 96-1
    lam = 0.001
    Map = np.ndarray([cn, 25, 26])
    Res = np.ndarray([cn, 25])
    meanPSNR = np.zeros(clusterszA)
    for t in tqdm(range(0, cn, 1)):
        c1 = Center[..., t].reshape(1, -1)
        # This part takes a most of the time so the RANSAC shouldn't add too much to that
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        # DO RANSAC HERE----------------------------------------------------------------------------------------
        # ------------------- RANSAC PARAMETERS -----------------------------------------------------------------
        indices = idx[0:clusterszA]
        numTrials = 200  # This bit takes 0.6% of the total time so 200 trials wont severely affect training time
        nRansac = 10  # Size of subset to take from indices for each iteration of ransac
        thRansac = 1.1  # Threshold for detecting outliers(PSNR_best/PSNR_current > thRansac then outlier)
        # -------------------------------------------------------------------------------------------------------
        # Output ransac variables
        minoutliers = np.inf
        Mbest = np.zeros([len(Center), len(Center)])
        PSNRbest = np.zeros(clusterszA)
        # Reference variables
        # Can just store means rather than storing LR1, HR1 and recalculating means
        LR1 = np.squeeze(X1[..., indices], axis=2)
        HR1 = np.squeeze(Y1[..., indices], axis=2)
        meanLR = LR1.mean(axis=0)
        meanHR = HR1.mean(axis=0)
        meanDiff = np.abs(meanHR-meanLR)
        LR1 = np.squeeze(X[..., indices], axis=2)
        LR1 = np.append(LR1, np.ones([1,LR1.shape[1]]), axis=0)
        for trial in range(0, numTrials):
            subsample = random.sample(range(0, len(indices)), nRansac)
            ind = indices[subsample]
            HR = np.squeeze(Y[..., ind], axis=2)
            LR = np.squeeze(X[..., ind], axis=2)
            LR = np.append(LR, np.ones([1,LR.shape[1]]), axis=0)
            M = HR.dot(LR.transpose().dot(inv(LR.dot(LR.transpose())
                                              + lam * np.identity(len(LR)))))

            SR = np.dot(M, LR1)
            SR1 = SR + np.tile(meanLR, [SR.shape[0], 1])

            mse = ((HR1 - SR1) ** 2).mean(axis=0)
            psnr = 10 * np.log10(mse ** -1)
            maxpsnr = psnr.max()
            noutliers = sum((maxpsnr / i) > thRansac for i in psnr)
            if (noutliers < minoutliers):
                minoutliers = noutliers
                Mbest = M
                PSNRbest = psnr
        # print("Min Outliers = {}".format(minoutliers))
        # fig1 = plt.figure()
        # fig1.suptitle("PSNR")
        # ax1 = fig1.add_subplot(111)
        # ax1.plot(np.arange(0,len(PSNRbest)), PSNRbest)
        # fig2 = plt.figure()
        # fig2.suptitle("Offsets")
        # ax2 = fig2.add_subplot(111)
        # ax2.plot(meanLR, linestyle=':')
        # ax2.plot(meanHR, linestyle=':')
        # ax2.plot(meanDiff, linestyle=':')
        # plt.show()
        meanPSNR = meanPSNR + PSNRbest
        # END RANSAC HERE---------------------------------------------------------------------------------------
        Map[t, :, :] = Mbest
        LR = np.squeeze(X[..., indices], axis=2)
        HR = np.squeeze(Y[..., indices], axis=2)
        # Res[t, :] = HR.mean(axis=1) - np.dot(Mbest, LR.mean(axis=1))
    meanPSNR = meanPSNR/cn
    fig = plt.figure()
    fig.suptitle("Mean PSNR across all runs")
    ax = fig.add_subplot(111)
    ax.plot(meanPSNR)
    plt.show()
    # -------------- SAVE MAPPINGS TO FILE ---------------------------------------------------------------------
    savfilename = 'data_files/pyMap' + str(i) + 'mat' + str(clusterszA + 1) + '.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return


def mapping_calculation_NoBias(dxloc, dyloc, i, clusterszA):
    xdict = scipy.io.loadmat(dxloc)
    ydict = scipy.io.loadmat(dyloc)
    X = xdict['X']
    Y = ydict['Y']
    loadfilename = 'data_files/pyCenter' + str(i) + '.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center']
    cn = len(Center[0])
    # clusterszA = 96-1
    lam = 0.01
    Map = np.ndarray([cn, 25, 25])
    Res = np.ndarray([cn, 25])
    for t in tqdm(range(0, cn, 1)):
        c1 = Center[..., t].reshape(1, -1)
        # This cdist part takes a very long time
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        idx = np.argsort(D, axis=0)
        LR = np.squeeze(X[..., idx[0:clusterszA]], axis=2)
        HR = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
        # # Add in bias term so regression model learns bias
        # LR = np.append(LR, np.ones([1,LR.shape[1]]), axis=0)
        M = HR.dot(LR.transpose().dot(inv(LR.dot(LR.transpose())
                                          + lam * np.identity(len(Center)))))
        Map[t, :, :] = M
        Res[t, :] = HR.mean(axis=1) - np.dot(M, LR.mean(axis=1))
    savfilename = 'data_files/pyMap' + str(i) + 'mat' + str(clusterszA + 1) + '.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return
