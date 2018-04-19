from run_kmeans import run_kmeans
from mapping_calculation import mapping_calculation
from cent_to_heir import cent_to_heir

dx = 'DX_and_DY/DX_all.mat'
dy = 'DX_and_DY/DY_all.mat'
run_kmeans(dx, 4096)
mapping_calculation(dx,dy,4096,96-1)
cent_to_heir(4096)