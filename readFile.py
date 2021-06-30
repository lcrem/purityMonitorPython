###################################################################################
#  Python script to read hdf5 from ProtoDUNE purity monitor runs
#  Author: Dr L. Cremonesi
#  Questions to l.cremonesi@qmul.ac.uk
###################################################################################

import h5py
import sys
import numpy as np

if len( sys.argv ) != 4:
    print('USAGE : %s <input folder > < run number > < PrM >')
    sys.exit (1)

inFolderName = sys.argv [1]
runNumber = int(sys.argv [2])
prm = int(sys.argv[3])

print( " Reading from " + inFolderName + "run " + str(runNumber) + " PrM " + str(prm))

f = h5py.File(inFolderName+'/Run'+str(runNumber)+'/PrM'+str(prm)+'filtered.hdf5', 'r')

groups=list(f.keys())

print(groups)

group = f.get('channel1')
pfitK = group.get('pfit')
perrK = group.get('perr')
group = f.get('channel2')
pfitA = group.get('pfit')
perrA = group.get('perr')


print(pfitK[:], pfitA[:])

QA=pfitA[3]
QK=pfitK[3]
R= -QA/QK
t1=pfitK[4]-pfitK[2]
t2=pfitA[2]-pfitK[4]
t3=pfitA[4]-pfitA[2]
taulife=-1/np.log(R)*(t2+0.5*(t1+t3))

print(' QA = %.1f mV \n QK = %.1f mV \n R = %.3f \n t1 = %.1f us \n t2 = %.1f us \n t3 = %.1f us \n taulife = %.1f us' % (QA, QK, R, t1*1e6, t2*1e6, t3*1e6, taulife*1e6))
