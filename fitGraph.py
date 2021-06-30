###################################################################################
#  Python script to fit RawAverages rootfiles from ProtoDUNE purity monitor runs
#  Author: Dr L. Cremonesi
#  Questions to l.cremonesi@qmul.ac.uk
###################################################################################

import sys, numbers, array, ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate
from scipy.optimize import curve_fit, leastsq
import h5py

if len( sys.argv ) != 5:
    print('USAGE : %s <input folder > < run number > < HV off run> < PrM >')
    sys.exit (1)

inFolderName = sys.argv [1]
runNumber = int(sys.argv [2])
hvoffNumber = int(sys.argv [3])
prm = int(sys.argv[4])

print( " Reading from " + inFolderName + "run " + str(runNumber) + " with HV off run " + str(hvoffNumber) + " PrM " + str(prm))

def getArraysFromFile(inFolderName, runNumber, channel, graphname):

    inFile1 = ROOT.TFile.Open ( inFolderName+'/Run'+str(runNumber)+'/RawAverages_ch'+str(channel)+'.root' ," READ ")
    graph1 = inFile1.Get(graphname)
    x = (np.array ( graph1.GetX() ) - 0.05*graph1.GetN()*2)*1e-9 
    y = np.array ( graph1.GetY() ) 

    # zero-baseline
    window_size=int(0.045*graph1.GetN())
    baseline = np.average(y[100:window_size])
    y = y - baseline
    
    return x,y


def greenFunction( t, F_elec, T_life, t0, Q0, t1, baseline):
    
    start=1./((t1-t0)*((1./T_life)-(1./F_elec)));
    fallExp=np.exp(-1*(t-t0)/F_elec);
    riseExp=np.exp(-1*(t-t0)/T_life);
    bigExp=fallExp*(1 - np.exp(-1*(t1-t0)*(1./T_life - 1./F_elec)));
    
    y = np.where(t<t0,  baseline, (np.where( t<t1, baseline+Q0*start*(fallExp-riseExp), baseline+Q0*start*bigExp ) ) )
    
    return y


if prm==1:
    channels=[1,2]
else:
    channels=[4,5]

fs = 5e8 # Sampling frequency
    
print ("Sampling frequency is " + str(fs) )
    
freq=200      # kHz
fc = freq*1000  # Cut-off frequency of the filter
print('Low-pass Butterworth filter ' + str(fc) + 'kHz')
w = fc / (fs / 2) # Normalize the frequency
b, a = signal.butter(4, w, 'low')

fout = h5py.File(inFolderName+'/Run'+str(runNumber)+'/PrM'+str(prm)+'filtered.hdf5', "w") 



for channel in channels:
    x1,y1 = getArraysFromFile(inFolderName, runNumber, channel, "justAvg")
    x1hvoff,y1hvoff = getArraysFromFile(inFolderName, hvoffNumber, channel, "justAvg")
    
    plt.figure(channel) 
    plt.plot(x1, y1, label='pre-filter')
    
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude [mV]")
    
    y1_butfilt = signal.filtfilt(b, a, y1);
    plt.plot(x1, y1_butfilt, label='lowpass '+str(freq)+'kHz')
    
    y1sub = np.subtract(y1,y1hvoff)
    #plt.plot(x1, y1sub, label='subtracted')
    y1_sub_filt = signal.filtfilt(b, a, y1sub);
    plt.plot(x1, y1_sub_filt, label='sub + lowpass '+str(freq)+'kHz')
    
    f = interpolate.interp1d(x1, y1_sub_filt, kind="linear")
    nbins = (x1[len(x1)-1]-x1[0])/0.1e-6 # interpolate every 0.1 us
    xnew = np.linspace(x1[0], x1[len(x1)-1], int(nbins)) 
    y_int = f(xnew)
    window_size, poly_order = 201, 3
    yy_sg = signal.savgol_filter(y_int, window_size, poly_order)
    plt.plot(xnew, yy_sg, label='sub +  Savitzky-Golay filter')    
    
    rms = np.sqrt(np.mean(y1hvoff[500000:]**2))
    print(rms)
    err_stdev = [rms for item in xnew]
    
    
    
    def fit_curvefit(p0, datax, datay, function, yerr=err_stdev, **kwargs):

        pfit, pcov = curve_fit(function,datax,datay,p0=p0,\
                                sigma=yerr, epsfcn=0.0001, **kwargs)
        error = [] 
        for i in range(len(pfit)):
            try:
                error.append(np.absolute(pcov[i][i])**0.5)
            except:
                error.append( 0.00 )
        pfit_curvefit = pfit
        perr_curvefit = np.array(error)
        return pfit_curvefit, perr_curvefit
    
    
    # fit curve
    if  channel==1 or channel==4 :
        Q0=np.amin(yy_sg)
        t0=0.6e-6
        t1=xnew[np.argmin(yy_sg)]
    else :
        Q0=np.max(yy_sg)
        t1=xnew[np.argmax(yy_sg)]
        t0=t1-10.e-6

    prefit = 250.e-6, 0.001, t0, Q0, t1, 1
    print("prefit = ", prefit)
    pfit, perr = fit_curvefit( prefit,  xnew, yy_sg, greenFunction, absolute_sigma=True)
    print("pfit = ", pfit)
    print("perr = ", perr)
    F_elec, T_life, t0, Q0, t1, baseline = pfit
    y_new = greenFunction(xnew, F_elec, T_life, t0, Q0, t1, baseline)
    print(' F_elec = %.5f us \n T_life = %.5f us \n t0 = %.5f us \n Q0 = %.5f \n t1 = %.5f us \n baseline = %.5f mV' % (F_elec*1e6, T_life*1e6, t0*1e6, Q0, t1*1e6, baseline))
    plt.plot(xnew, y_new, label='fitted')

    grp = fout.create_group('channel'+str(channel))
    grp.create_dataset('x_int', data=xnew)
    grp.create_dataset('y_int', data=y_new)
    grp.create_dataset('pfit',  data=pfit)
    grp.create_dataset('perr',  data=perr)

    plt.legend()

plt.show()
fout.close()
