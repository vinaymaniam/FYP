load data_files/pyCenter16384.mat

tic
[heirarchy,index] = c2h(Center);
toc

save data_files/newheir16384.mat heirarchy index