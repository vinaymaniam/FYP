import scipy.io
from sklearn.cluster import MiniBatchKMeans
import time

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
    return
