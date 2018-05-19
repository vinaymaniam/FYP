from run_kmeans import run_kmeans
from mapping_calculation import mapping_calculation, mapping_calculationRANSAC
from cent_to_heir import cent_to_heir
import matlab.engine

def runFullTraining(rkm=1,rmc=1,rcm2c=1):
    dx = 'DX_and_DY/DX_all.mat'
    dy = 'DX_and_DY/DY_all.mat'

    n = 4096
    success = 1
    if(rkm==1):
        run_kmeans(dx, n)
    if(rmc==1):
        # mapping_calculation(dx,dy,n,96-1)
        mapping_calculationRANSAC(dx, dy, n, 96 - 1)
    if(rcm2c==1):
        cent_to_heir(n)
        eng = matlab.engine.start_matlab()
        success = eng.ConvMat2Cell(n)
    return success

runFullTraining(0,1,1)