import numpy as np
import time

a = np.random.rand(25,10000000)
b = np.random.rand(25,25)

l = np.int32(np.linspace(100,10000000,20))
tvec = np.zeros([20,1])
iter = 0
for i in l:
    a1 = a[:,0:i]
    t1 = time.time()
    c = b.dot(a1)
    t2 = time.time()
    tvec[iter] = t2-t1
    print(tvec[iter])
    iter = iter+1

print('done')