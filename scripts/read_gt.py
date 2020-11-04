import numpy as np
def get_gt(filename,shape):
    dt=np.loadtxt(filename)
    dt=dt[:,[3,7,11]]
    dt=np.expand_dims(dt,0)-np.expand_dims(dt,1)
    dt=np.sqrt(np.sum(np.square(dt),axis=-1))
    m=dt.shape[0]//shape
    dt=dt[::m,::m]
    dt=dt[:shape,:shape]
    return dt

