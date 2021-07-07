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


errQA=perrA[3]
errQK=perrK[3]
errR= R*np.sqrt( (errQA/QA)*(errQA/QA) + (errQK/QK)*(errQK/QK) )
errt1=np.sqrt(perrK[4]*perrK[4] + perrK[2]*perrK[2])
errt2=np.sqrt(perrA[2]*perrA[2] + perrK[4]*perrK[4])
errt3=np.sqrt(perrA[4]*perrA[4] + perrA[2]*perrA[2])
errtaulife=0


print(' QA =  %3.1f +/- %3.1f mV \n QK = %3.1f +/- %3.1f mV \n R  =  %3.2f +/- %3.2f \n t1 =  %3.1f +/- %3.1f us \n t2 =  %3.1f +/- %3.1f us \n t3 =  %3.1f +/- %3.1f us \n taulife =  %3.1f +/- %3.1f us' % (QA, errQA, QK, errQK, R, errR, t1*1e6, errt1*1e6, t2*1e6, errt2*1e6, t3*1e6, errt3*1e6, taulife*1e6, errtaulife*1e6))
