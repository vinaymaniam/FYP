import numpy as np
from numpy.linalg import inv
import scipy
import scipy.io
import time

xdict = scipy.io.loadmat('DX_all.mat')
ydict = scipy.io.loadmat('DY_all.mat')
X = xdict['X']
Y = ydict['Y']
i=8192
# -----------------------------------
tstart = time.time()
loadfilename = 'pyCenter'+str(i)+'.mat'
cent = scipy.io.loadmat(loadfilename)
Center = cent['Center']
cn = len(Center[0])
clusterszA = 96-1
lam = 0.01
Map = np.ndarray([cn,25,25])
for t in range(0,cn,1):
    c1 = Center[...,t].reshape(-1,1)
    D = scipy.spatial.distance.cdist(X.transpose(), Center[...,t].reshape(1,-1))
    # This sorting part takes a long time!
    dv = D.sort()
    idx = D.argsort()
    patchesL = np.squeeze(X[...,idx[0:clusterszA]], axis=2)
    patchesH = np.squeeze(Y[..., idx[0:clusterszA]], axis=2)
    M = patchesH.dot(patchesL.transpose().dot(
        inv(patchesL.dot(patchesL.transpose())
            + lam*np.identity(len(Center)))))
    Map[t,:,:] = M
    if(np.mod(t,32)==0):
        print('Finished ',t,' out of ',cn)
tend = time.time()
print('Map generation took ',(tend-tstart),' seconds')
savfilename = 'pyMap'+str(i)+'cell'+str(clusterszA+1)+'.mat'
data = {'Center': Map}
scipy.io.savemat(savfilename, data)
