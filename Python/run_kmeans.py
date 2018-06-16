import scipy.io
from sklearn.cluster import MiniBatchKMeans, KMeans
import numpy as np
import time
import matplotlib.pyplot as plt

def run_kmeans(loadfilename, nctrds):
    Xdict = scipy.io.loadmat(loadfilename)
    X = Xdict['X'].transpose()
    t1 = time.time()
    # kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', init_size=2*nctrds, max_iter=500, verbose=False).fit(X)
    kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', batch_size=10000, max_iter=500, verbose=0, n_init=3).fit(X)
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
    return

def run_kmeans_nonminibatch(loadfilename, nctrds):
    Xdict = scipy.io.loadmat(loadfilename)
    X = Xdict['X'].transpose()
    print('About to run K-means for n= '+str(nctrds))
    t1 = time.time()
    kmeans = KMeans(n_clusters=nctrds, init='k-means++', max_iter=500, verbose=1, n_init=1).fit(X)
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
    return
