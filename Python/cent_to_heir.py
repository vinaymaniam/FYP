import numpy as np
import scipy
import scipy.io
import time
import sklearn
from sklearn.cluster import MiniBatchKMeans

sz = 8192
loadfilename = 'pyCenter'+str(sz)+'.mat'
cent = scipy.io.loadmat('pyCenter8192.mat')
Center = cent['Center'].transpose()
nctrds = int(round(np.sqrt(len(Center))))

kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', init_size=2 * nctrds, max_iter=500, verbose=False).fit(
    Center)
C = kmeans.cluster_centers_;

knn = sklearn.neighbors.NearestNeighbors(n_neighbors=3 * nctrds, metric='sqeuclidean')
knn.fit(Center)
index = knn.kneighbors(C)[1]
heirarchy = np.ndarray([nctrds, 25, len(index[0]) + 1])
heirarchy[:, :, 0] = C
for i in range(0, len(index[0]), 1):
    for j in range(0, nctrds, 1):
        heirarchy[j, :, i + 1] = Center[index[j, i], :]
data = {'heirarchy': heirarchy, 'index': index}
savfilename = 'pyHeirarchy'+str(sz)+'.mat'
scipy.io.savemat(savfilename, data)
