import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import numpy as np
from read_gt import *
import sys
z= np.loadtxt('OUT.txt')
gt=get_gt('/mnt/external/kitti-devkit-odom/ground_truth/poses/'+sys.argv[1]+'.txt',z.shape[0])
for i in range(z.shape[0]):
    z[i,i+1:]=1.0
    z[i,i]=0.0
    ind = np.nonzero(gt[i,i+5:]<4)
    z[i,i+5:][ind]=0.0
z=-ndimage.morphology.grey_dilation(-z,size=(3,3))
clip=0.04
z[z>clip]=clip
plt.imshow(z,cmap=plt.get_cmap('bone'))
plt.xlabel("Pointcloud index")
plt.ylabel("Pointcloud index")
plt.title("Similarity matrix.")
plt.legend()
plt.show()
