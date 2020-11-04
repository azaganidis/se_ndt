import numpy as np
dt=np.fromfile('loop.dat',dtype=np.float32)
num=int(np.sqrt(dt.shape[0]*2))
z=np.zeros((num,num))
start=0
for i in range(num-1):
    start=start+i
    z[i,:i+1]=dt[start:start+i+1]
np.savetxt('OUT.txt',z)
