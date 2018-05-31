import scipy.io
from sklearn.cluster import MiniBatchKMeans
import numpy as np
import time
import matplotlib.pyplot as plt

def run_kmeans(loadfilename, nctrds):
    Xdict = scipy.io.loadmat(loadfilename)
    X = Xdict['X'].transpose()
    t1 = time.time()
    # nctrds = 4096
    kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', init_size=2*nctrds, max_iter=500, verbose=False).fit(X)
    t2 = time.time()
    data = {'Center': kmeans.cluster_centers_.transpose()}

    savefilename = 'data_files/pyCenter' + str(nctrds) + '.mat'
    scipy.io.savemat(savefilename,data)
    print('K-means took ',(t2-t1),' seconds')
    print('Inertia = ',kmeans.inertia_)
    # Compute how many of elements in each cluster
    clustercount = np.zeros(nctrds)
    for i in range(0,len(X)):
        clustercount[kmeans.labels_[i]] = clustercount[kmeans.labels_[i]] + 1
    print('Number of centroids with 0 or 1 elements = ' + str(np.count_nonzero(clustercount==0)+np.count_nonzero(clustercount==1)))
    # fig1 = plt.figure()
    # fig1.suptitle("Cluster Count For Each Centroid")
    # ax1 = fig1.add_subplot(111)
    # ax1.plot(np.arange(0, len(clustercount)), clustercount)
    # plt.show()
    return
