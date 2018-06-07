from run_kmeans import run_kmeans
# from mapping_calculation import mapping_calculation, mapping_calculationRANSAC
from mapping_calc_pre_crossval import mapping_calculation, mapping_calculation_3by3, mapping_calculation_NbyN
from cent_to_heir import cent_to_heir, cent_to_heir_NbyN
import matlab.engine
import numpy as np

def runFullTraining(rkm=1,rmc=1,rcm2c=1):
    stage = 4
    n = 16384
    # numneighbors = 96 - 1
    # nn = np.array([12,48,96,192,384,768,3072]) - 1
    nn = np.array([192]) - 1
    success = 1
    mode = 0 # 0:5x5, 3:3x3, 4:4x4 etc
    # for n in [128,256,512,1024,2048,4096,8192,16384]:
    for n in [2048]:
        if(mode == 0):
            for numneighbors in nn :
                if(stage == 1):
                    dx = 'DX_and_DY/DX_all.mat'
                    dy = 'DX_and_DY/DY_all.mat'
                elif (stage == 2):
                    dx = 'DX_and_DY/DX_all'+str(n)+'.mat'
                    dy = 'DX_and_DY/DY_all'+str(n)+'.mat'
                elif (stage >= 3):
                    dx = 'DX_and_DY/DX_all' + str(n) + '_' + str(stage-1) + '.mat'
                    dy = 'DX_and_DY/DY_all' + str(n) + '_' + str(stage-1) + '.mat'

                if(rkm==1):
                    run_kmeans(dx, n)
                    cent_to_heir(n, stage)
                if(rmc==1):
                    # ADDED BIAS TERM TO mapping_calculation(not yet to ransac)
                    # CHANGE: MADE clusterszA inversely proportional to nCtrds
                    # numneighbors = int(np.rint(233013/n))
                    mapping_calculation(dx,dy,n,numneighbors)
                    # mapping_calculationRANSAC(dx, dy, n, 96 - 1)
                if(rcm2c==1):
                    eng = matlab.engine.start_matlab()
                    success = eng.ConvMat2Cell(n, stage, int(numneighbors+1))
                print('Finished training model for clusterszA = ' + str(numneighbors))
        else:
            numneighbors = int(nn)
            if (stage == 1):
                dx = 'DX_and_DY/DX_all{}x{}.mat'.format(mode,mode)
                dy = 'DX_and_DY/DY_all{}x{}.mat'.format(mode,mode)

            if (rkm == 1):
                run_kmeans(dx, n)
                cent_to_heir_NbyN(n, stage, mode)
            if (rmc == 1):
                # ADDED BIAS TERM TO mapping_calculation(not yet to ransac)
                # CHANGE: MADE clusterszA inversely proportional to nCtrds
                # numneighbors = int(np.rint(233013/n))
                mapping_calculation_NbyN(dx, dy, n, numneighbors,mode)
                # mapping_calculationRANSAC(dx, dy, n, 96 - 1)
            if (rcm2c == 1):
                eng = matlab.engine.start_matlab()
                success = eng.ConvMat2CellNxN(n, stage, int(numneighbors + 1), mode)
            print('Finished training model for clusterszA = ' + str(numneighbors))
    return success

runFullTraining(1,1,1)