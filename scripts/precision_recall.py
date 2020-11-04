from sklearn.metrics import roc_curve, precision_recall_curve
from matplotlib import pyplot
import numpy as np
from read_gt import *
def pr(dt,z,dist_thr):
    skip=5
    GT=[]
    EST=[]
    for i in range(skip,z.shape[0]):
        gt = dt[i,:i-skip+1]
        st = z[i,:i-skip+1]
        GT+=(gt<dist_thr).tolist()
        EST+=st.tolist()
    GT=np.array(GT)
    EST=np.array(EST)
    fpr,tpr,_=roc_curve(GT,EST)
    pyplot.plot(fpr,tpr,marker='.',label='SE-NDT')
    pyplot.xlabel('False positive rate')
    pyplot.ylabel('True positive rate')
    pyplot.legend()
    pyplot.show()
    precision,recall,_=precision_recall_curve(GT,EST)
    pyplot.plot(recall,precision,marker='.',label='SE-NDT')
    pyplot.xlabel('Recall')
    pyplot.ylabel('Precision')
    pyplot.legend()
    pyplot.show()
z=np.loadtxt('OUT.txt')
dt=get_gt('/mnt/external/kitti-devkit-odom/ground_truth/poses/00.txt',z.shape[0])
pr(dt,z,10)

