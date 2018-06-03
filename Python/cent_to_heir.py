import numpy as np
import scipy
import scipy.io
import time
import sklearn
from sklearn.cluster import MiniBatchKMeans


def cent_to_heir(sz, stage):
    loadfilename = 'data_files/pyCenter'+str(sz)+'.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center'].transpose()
    nctrds = int(round(np.sqrt(len(Center))))

    kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', init_size=2 * nctrds, max_iter=500, verbose=False).fit(
        Center)
    C = kmeans.cluster_centers_

    knn = sklearn.neighbors.NearestNeighbors(n_neighbors=3 * nctrds, metric='sqeuclidean')
    knn.fit(Center)
    index = knn.kneighbors(C)[1]
    heirarchy = np.ndarray([nctrds, 25, len(index[0]) + 1])
    heirarchy[:, :, 0] = C
    for i in range(0, len(index[0]), 1):
        for j in range(0, nctrds, 1):
            heirarchy[j, :, i + 1] = Center[index[j, i], :]

    index = index + 1 # matlab requires 1 indexed arrays
    data = {'heirarchy': heirarchy, 'index': index}
    savfilename = 'data_files/'+str(stage)+'pyHeirarchy'+str(sz)+'.mat'
    scipy.io.savemat(savfilename, data)

def cent_to_heir_NbyN(sz, stage, mode):
    loadfilename = 'data_files/pyCenter'+str(sz)+'.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center'].transpose()
    nctrds = int(round(np.sqrt(len(Center))))

    kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', init_size=2 * nctrds, max_iter=500, verbose=False).fit(
        Center)
    C = kmeans.cluster_centers_

    knn = sklearn.neighbors.NearestNeighbors(n_neighbors=3 * nctrds, metric='sqeuclidean')
    knn.fit(Center)
    index = knn.kneighbors(C)[1]
    heirarchy = np.ndarray([nctrds, mode*mode, len(index[0]) + 1])
    heirarchy[:, :, 0] = C
    for i in range(0, len(index[0]), 1):
        for j in range(0, nctrds, 1):
            heirarchy[j, :, i + 1] = Center[index[j, i], :]

    index = index + 1 # matlab requires 1 indexed arrays
    data = {'heirarchy': heirarchy, 'index': index}
    savfilename = 'data_files/'+str(stage)+'pyHeirarchy'+str(sz)+'_{}x{}.mat'.format(mode,mode)
    scipy.io.savemat(savfilename, data)




def cent_to_heir_3by3(sz, stage):
    loadfilename = 'data_files/pyCenter'+str(sz)+'.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center'].transpose()
    nctrds = int(round(np.sqrt(len(Center))))

    kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', init_size=2 * nctrds, max_iter=500, verbose=False).fit(
        Center)
    C = kmeans.cluster_centers_

    knn = sklearn.neighbors.NearestNeighbors(n_neighbors=3 * nctrds, metric='sqeuclidean')
    knn.fit(Center)
    index = knn.kneighbors(C)[1]
    heirarchy = np.ndarray([nctrds, 9, len(index[0]) + 1])
    heirarchy[:, :, 0] = C
    for i in range(0, len(index[0]), 1):
        for j in range(0, nctrds, 1):
            heirarchy[j, :, i + 1] = Center[index[j, i], :]

    index = index + 1 # matlab requires 1 indexed arrays
    data = {'heirarchy': heirarchy, 'index': index}
    savfilename = 'data_files/'+str(stage)+'pyHeirarchy'+str(sz)+'_3x3.mat'
    scipy.io.savemat(savfilename, data)

def cent_to_heir_6by6(sz, stage):
    loadfilename = 'data_files/pyCenter'+str(sz)+'.mat'
    cent = scipy.io.loadmat(loadfilename)
    Center = cent['Center'].transpose()
    nctrds = int(round(np.sqrt(len(Center))))

    kmeans = MiniBatchKMeans(n_clusters=nctrds, init='k-means++', init_size=2 * nctrds, max_iter=500, verbose=False).fit(
        Center)
    C = kmeans.cluster_centers_

    knn = sklearn.neighbors.NearestNeighbors(n_neighbors=3 * nctrds, metric='sqeuclidean')
    knn.fit(Center)
    index = knn.kneighbors(C)[1]
    heirarchy = np.ndarray([nctrds, 36, len(index[0]) + 1])
    heirarchy[:, :, 0] = C
    for i in range(0, len(index[0]), 1):
        for j in range(0, nctrds, 1):
            heirarchy[j, :, i + 1] = Center[index[j, i], :]

    index = index + 1 # matlab requires 1 indexed arrays
    data = {'heirarchy': heirarchy, 'index': index}
    savfilename = 'data_files/'+str(stage)+'pyHeirarchy'+str(sz)+'_6x6.mat'
    scipy.io.savemat(savfilename, data)