from run_kmeans import run_kmeans
# from mapping_calculation import mapping_calculation, mapping_calculationRANSAC
from mapping_calc_pre_crossval import mapping_calculation, mapping_calculation_quad
from cent_to_heir import cent_to_heir
import matlab.engine
import numpy as np

def runFullTraining(rkm=1,rmc=1,rcm2c=1):
    stage = 1
    n = 16384
    success = 1
    if(stage == 1):
        dx = 'DX_and_DY/DX_all.mat'
        dy = 'DX_and_DY/DY_all.mat'
    else:
        dx = 'DX_and_DY/DX_all'+str(n)+'.mat'
        dy = 'DX_and_DY/DY_all'+str(n)+'.mat'

    if(rkm==1):
        run_kmeans(dx, n)
    if(rmc==1):
         # ADDED BIAS TERM TO mapping_calculation(not yet to ransac)
        # CHANGE: MADE clusterszA inversely proportional to nCtrds
        numneighbors = 96-1
        # numneighbors = int(np.rint(233013/n))
        mapping_calculation(dx,dy,n,numneighbors)
        # mapping_calculationRANSAC(dx, dy, n, 96 - 1)
    if(rcm2c==1):
        cent_to_heir(n, stage)
        eng = matlab.engine.start_matlab()
        success = eng.ConvMat2Cell(n, stage)
    return success

runFullTraining(1,1,1)