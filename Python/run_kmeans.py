import scipy.io
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
import numpy as np
import time

Xdict = scipy.io.loadmat('DX_all.mat')
X = Xdict['X'].transpose()
t1 = time.time()
# kmeans = KMeans(n_clusters=4, init='k-means++', precompute_distances=False).fit(X)
numcentroids = 4096
kmeans = MiniBatchKMeans(n_clusters=numcentroids, init='k-means++', init_size=2*numcentroids, max_iter=250, verbose=False).fit(X)
t2 = time.time()
data = {'Center': kmeans.cluster_centers_}
scipy.io.savemat('pyCenter4096.mat',data)
print('K-means took ',(t2-t1),' seconds')
print('Inertia = ',kmeans.inertia_)
