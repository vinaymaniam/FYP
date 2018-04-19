import numpy as np
from numpy.linalg import inv
import scipy
import scipy.io
import time
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
        D = scipy.spatial.distance.cdist(X.transpose(), c1)
        # This sorting part takes a long time!
        idx = np.argsort(D, axis=0)
        LR = np.squeeze(X[..., idx[0:clusterszA]], axis=2)
        HR = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
        M = HR.dot(LR.transpose().dot(inv(LR.dot(LR.transpose())
                   + lam * np.identity(len(Center)))))
        Map[t,:,:] = M
        Res[t,:] = HR.mean(axis=1) - np.dot(M,LR.mean(axis=1))
    savfilename = 'data_files/pyMap'+str(i)+'cell'+str(clusterszA+1)+'.mat'
    data = {'Map': Map, 'Res': Res}
    scipy.io.savemat(savfilename, data)
    return
