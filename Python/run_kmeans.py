import scipy.io
from sklearn.cluster import MiniBatchKMeans
import time

Xdict = scipy.io.loadmat('DX_all.mat')
X = Xdict['X'].transpose()
t1 = time.time()
numcentroids = 32768
kmeans = MiniBatchKMeans(n_clusters=numcentroids, init='k-means++', init_size=2*numcentroids, max_iter=500, verbose=False).fit(X)
t2 = time.time()
data = {'Center': kmeans.cluster_centers_.transpose()}

savfilename = 'pyCenter' + str(numcentroids) + '.mat'
scipy.io.savemat(savfilename,data)
print('K-means took ',(t2-t1),' seconds')
print('Inertia = ',kmeans.inertia_)