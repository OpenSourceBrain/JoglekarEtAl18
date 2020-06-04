#
# Work in progress ...
#
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib
from scipy import signal

# To do
# Peak gamma. You should change the input on V1

def parameters(areas):
    
    
    nAreas=len(areas)
    dt=0.2e-3
    triallength=20
    transient=4
    sc2=1
    sc5=5
    tau=np.array([0.006*sc2, 0.015*sc2, 0.006*sc5, 0.015*sc5])
    tstep=(np.ones((4,nAreas)).T*(dt/tau)).T
    sig=0*np.array([0.3, 0.3, 0.45, 0.45]) # there is no noise
    tstep2=((dt*sig*sig)/tau)**0.5
    tstep2=(np.ones((4,nAreas)).T*tstep2).T
    binx=20
    eta=0.2
    G=1.1
    
    # local and interlaminar coupling
    
    J2e=1.00   # L2/3 excit to L5 excit,
    J2i=0.     # L2/3 excit to L5 inhib
    J5e=0.     # L5 excit to L2/3 excit
    J5i=0.75   # L5 excit to L2/3 inhib,

    J=np.zeros((4,4));
    
    #local, layer 2:
    
    J[0,0]=1.5;
    J[0,1]=-3.25;
    J[1,0]=3.5;
    J[1,1]=-2.5;

    #local, layer 5:
    J[2,2]=1.5;
    J[2,3]=-3.25;
    J[3,2]=3.5;
    J[3,3]=-2.5;

    #inter-laminar:
    J[2,0]=J2e;
    J[3,0]=J2i;
    J[0,2]=J5e;
    J[1,2]=J5i;
    
    inputbg=np.zeros((4,nAreas));
    inputbg[[0, 2],:]=0; #thalamus!
    
    data=loadmat('subgraphData30.mat')
    fln=data['flnMat']
    sln=data['slnMat']
    fln=fln[0:nAreas,0:nAreas]
    sln=sln[0:nAreas,0:nAreas]
    
    # Compreession
    fln=1.2*fln**0.3
    #selectivity matrix
    s=0.1*np.ones((nAreas,nAreas))
    
    Wff=fln*sln
    Wfb=fln*(np.ones((nAreas,nAreas))-sln)
    
    wiring=loadmat('subgraphWiring30.mat')
    wires=wiring['wiring']
    wires=wires[0:nAreas,0:nAreas]
    delay= np.round(wires/(1500*dt)+1).astype(int)
    
    #we properly normalize the FLNs entering each area:
    normff=np.sum(Wff,axis=1);
    normfb=np.sum(Wfb,axis=1);
    for i in range(nAreas):
        Wff[i,:]=G*(Wff[i,:]/normff[i]);
        Wfb[i,:]=G*(Wfb[i,:]/normfb[i]);

    # Dictionary with parameters
    par={'dt':dt,
         'triallength':triallength,
         'transient':transient,
         'dt':dt,
         'tau':tau,
         'tstep':tstep,
         'tstep2':tstep2,
         'binx':binx,
         'eta':eta,
         'J':J,
         'inputbg':inputbg,
         's':s,
         'Wff':Wff,
         'Wfb':Wfb,
         'delay':delay
        }    

    return par


def hierarchy(indice,areas,par):

     nAreas=len(areas)
     
     nobs=len(np.arange(par['dt'],(par['triallength']-par['transient']),par['dt']))
     X=np.zeros((nAreas,round(nobs/par['binx'])))
     X2=np.zeros((nAreas,round(nobs/par['binx'])))
     X5=np.zeros((nAreas,round(nobs/par['binx'])))
     

     # external input on the full brain network:
     Iext=np.zeros((4,nAreas))
     Iext[0,0]=150 #input to V1


     # real simulation:
     if indice==0:   
        firingRate=trialdelays(0,par,Iext,nAreas)
     elif indice==1:
        firingRate=trialdelays(1,par,Iext,nAreas)

     t0=round((par['dt']+par['transient'])/par['dt'])
     t=np.arange(t0,max(np.shape(firingRate)),par['binx'])

     for i in range(nAreas):
         
        #we save the excit. firing rate re(L2/3)+re(L5) :
        X[i,:]=firingRate[0,t,i]*par['eta']+firingRate[2,t,i]*(1-par['eta'])
        X2[i,:]=firingRate[0,t,i]
        X5[i,:]=firingRate[2,t,i]
    
     return X,X2,X5


def trialdelays(iteration,par,Iext,nAreas):

    compls=np.ones((nAreas,nAreas))-par['s']
    auxones=np.ones((4,nAreas))

    Wfbe1=par['s']*par['Wfb']
    Wfbe2=compls*par['Wfb']
    
    J=par['J']
    
    if iteration==0:
        lrEtoE=1
        localItoE=1
    elif iteration==1:       
        lrEtoE=2.2
        localItoE=1.2
        J[0,1]=J[0,1]*localItoE;
 
    # Set variables
    totalinput=np.zeros((4,nAreas))
    irate=np.zeros((4,nAreas))
    iratenew=np.zeros((4,nAreas))
    inoise=np.zeros((4,nAreas))
    
    
    rate=np.zeros((4,round(par['triallength']/par['dt']),nAreas))
    transfer=np.zeros((4,nAreas))
    #noise for re2,ri2,re5,ri5
    xi=np.random.normal(0,1,(4,round(par['triallength']/par['dt']),nAreas)) 
    drate1=np.zeros((nAreas,nAreas))
    drate3=np.zeros((nAreas,nAreas))
    
    #we build the identity-block matrix:
    blockmatrix=np.zeros((nAreas*nAreas,nAreas))
    
    for i in range(nAreas):
        blockmatrix[(i*nAreas):((i+1)*nAreas),i]=1
    
    
    
    #first iteration:
    rate[:,0,:]=5*(1+np.tanh(2*xi[:,0,:])); # between 0 and 10 spikes/s

    #Now we start the real simulation:
    i=1    
    
    times=np.arange(par['dt'],par['triallength'],par['dt'])
    dprov2=np.zeros((nAreas*nAreas))
    for time in times:
 	
        #we set the instantaneous and delayed rates for computations:
        irate[:,:]=rate[:,i-1,:]
      
        #we set the matrix with the delays for this iteration:
        if time>0.2: #little transient 
            
            delaynow=(i-1)-par['delay']
            
            #for rates from L2/3e
            dprov=np.squeeze(rate[0,delaynow.flatten('F'),:])
            dprov=dprov*blockmatrix
            for idx in range(nAreas):
                dprov2[(idx*nAreas):((idx+1)*nAreas)]=dprov[(idx*nAreas):((idx+1)*nAreas),idx]
            drate1=np.reshape(dprov2,(nAreas,nAreas),order='F')      # mexer aqui
          
          
            #and for rates from L5e:
            dprov=np.squeeze(rate[2,delaynow.flatten('F'),:]);
            dprov=dprov*blockmatrix
            for idx in range(nAreas):
                dprov2[(idx*nAreas):((idx+1)*nAreas)]=dprov[(idx*nAreas):((idx+1)*nAreas),idx]
            drate3=np.reshape(dprov2,(nAreas,nAreas),order='F')      

     
        pulse=np.zeros((4,nAreas))
        if time>12 and time<12.3:
            pulse=Iext
  
    
        inputbg=par['inputbg']+pulse
        totalinput=inputbg+np.dot(J,irate)
  
        # interareal FF projections:
        totalinput[0,:]=totalinput[0,:]+1.0*lrEtoE*(np.sum(drate1*par['Wff'],axis=1)).T
        totalinput[0,:]=totalinput[0,:]+0.2*(np.sum(drate1*par['Wff'],axis=1)).T
        # interareal FB projections:
        totalinput[0,:]=totalinput[0,:]+1.0*(np.sum(drate3*Wfbe1,axis=1)).T
        totalinput[1,:]=totalinput[1,:]+0.5*(np.sum(drate3*par['Wfb'],axis=1)).T
        totalinput[2,:]=totalinput[2,:]+1.0*(np.sum(drate3*Wfbe2,axis=1)).T
        totalinput[3,:]=totalinput[3,:]+0.5*(np.sum(drate3*par['Wfb'],axis=1)).T
  
        # input after transfer functions:
        transfer=totalinput/(auxones-np.exp(-1*totalinput))

        # we update the firing rates of all areas:
        inoise=xi[:,i-1,:]
        iratenew=irate+par['tstep']*(-irate+transfer)+par['tstep2']*inoise
        
        rate[:,i,:]=iratenew

        #index iteration
        i=i+1;
    
    
    return rate


def figProp(X2,par,nAreas):
    
    ratios=np.zeros(nAreas);
    
    Tmin=7.5
    Tmax=9
   
    dt=par['binx']*par['dt']  
    
    for i in range(nAreas):
      
      r2=X2[i,round(Tmin/dt):round(Tmax/dt)]
      ratios[i]=r2[200]-r2[300]
      
    return ratios

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')


def peakGamma(X,fs,fmin):
    
    nAreas=np.shape(X)[0]
    
    peak=np.zeros(nAreas)
    
    for i in range(nAreas):
        
        f, Pxx_den = signal.periodogram(X[i,:], fs)
        pxx=movingaverage(Pxx_den, 400)
        z=np.where(f>fmin)
        
        newpxx=pxx[z]
        newf=f[z]
        
        pks,_=signal.find_peaks(newpxx)
        z=len(pks) # if there's at least one peak:
  
        if z>=1:
            z3=max(newpxx[pks])
            # frequency=newfx(loc(z3,1));  # location, in Hz, of the highest peak
            # amplitudeA=newpxx(loc(z3,1)); # power of the peak in the spectrum
        else:
            frequency=0;  # no oscillations
            amplitudeA=0;
      
        
    return peak


areas=np.arange(29)
par=parameters(areas)

X,X2,X5=hierarchy(0,areas,par)

rat0=figProp(X2,par,len(areas))

X,X2,X5=hierarchy(1,areas,par)
rat1=figProp(X2,par,len(areas))

plt.semilogy(areas,rat0,'green',marker=".",  markersize=5)
plt.semilogy(areas,rat1,'purple',marker=".",  markersize=5)
plt.show()

