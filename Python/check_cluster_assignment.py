import numpy as np
import scipy.io as spio
import scipy as sp
from tqdm import tqdm
import matplotlib.pyplot as plt

xdict = spio.loadmat('DX_and_DY/DX_all.mat')
X = xdict['X']
cent = spio.loadmat('data_files/pyCenter16384.mat')
center = cent['Center']

clustercount = np.zeros(len(center[0]))
for i in tqdm(range(0,len(X[0]))):
    x = X[:, i].reshape(1, -1)
    D = sp.spatial.distance.cdist(center.T,x)
    idx = np.argsort(D, axis=0)
    clustercount[idx[0]] = clustercount[idx[0]] + 1

fig1 = plt.figure()
fig1.suptitle("Cluster Count For Each Centroid")
ax1 = fig1.add_subplot(111)
ax1.plot(np.arange(0,len(clustercount)), clustercount)
plt.show()

print('done')