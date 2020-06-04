import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib

def threshLinear(x,limiar):
        return np.maximum(x,limiar)

def rateModel(ini,wEE,wEI,wIE,wII,wEElr,wIElr,nSteps,dt,noise,areaidx,areaInput,bkg,tauE,tauI):
    
    N=np.shape(wEElr)[0]
    bkg=np.squeeze(bkg)
    # initial conditions
    stateExc=ini[:N] 
    stateInh=ini[N:]
    
    # Store rates
    ratesExc=np.zeros((nNodes,int(nSteps)))
    ratesInh=np.zeros((nNodes,int(nSteps)))
    
    ratesExc[:,0]=stateExc
    ratesInh[:,0]=stateInh
    
    # Noise std can vary between areas.
    noiseStd=noise/np.sqrt(dt)
    
    # Simulation of rate model
    toy=0
    for i in range(1,int(nSteps)):
        
        if (areaInput[i]==0 and toy==0):
            stateExc = 10*np.ones(nNodes)
            stateInh = 35*np.ones(nNodes)
        else:    
            toy=1
        
            stateExcI=np.dot(wEElr,stateExc)+wEE*stateExc+ wEI*stateInh+bkg[:N]+noise*np.random.randn(nNodes)
            stateInhI=np.dot(wIElr,stateExc)+wIE*stateExc+ wII*stateInh+bkg[N:]+noise*np.random.randn(nNodes)
    
            # Add in the area specific input (single area, exc population)
            stateExcI[areaidx]=stateExcI[areaidx]+areaInput[i]
        
            stateExc = threshLinear(stateExc+(dt/tauE)*(-stateExc+threshLinear(stateExcI,10)),0);
            stateInh = threshLinear(stateInh+(dt/tauI)*(-stateInh+threshLinear(stateInhI,35)),0);
            
        ratesExc[:,i]=stateExc
        ratesInh[:,i]=stateInh;
    
    return ratesExc,ratesInh



# Parameters 

# anatomical
data=loadmat('subgraphData.mat')
nNodes=data['nNodes'][0,0]
hier=np.squeeze(data['hierVals'])
hierNorm=hier/max(hier)
fln=data['flnMat']

# Input 
tStart=2      
tEnd=2.25
ampl=21.8*1.9
nSteps=5e4
dt=2e-4
areaInput=np.zeros((int(nSteps)))
areaInput[round(tStart/dt) : round(tEnd/dt) ]=ampl
alpha=.68
betaE = .066
betaI = .351
omegaEE = 24.3
omegaEI = 19.7
omegaIE = 12.2
omegaII = 12.5
muIE = 25.3 
muEErange = np.arange(28,44,2)
noise=0
# Values or rate at steady state
desiredSs=np.concatenate((10*np.ones((nNodes,1)), 35*np.ones((nNodes,1)))) 

tauE=2e-2
tauI=1e-2

prop24c = np.zeros((len(muEErange),2));

areaidx=1


toplotstart=round(1.75/dt);
toplotend=round(5/dt);

# Store data for figure B
Area1data=np.zeros((2,int(nSteps)))
Area2data=np.zeros((2,int(nSteps)))
# Store data for figure C
Area2peak=np.zeros((2,len(muEErange)))


idxMu=0
for u in range(2):
    if u==1:
        fln=np.tril(fln)    
    for muIdx in range(len(muEErange)): #len(muEErange)
       
        muEE=muEErange[muIdx]
                  
        bgInh=np.zeros(nNodes)
        
        # Synaptic weights for intra-areal connections
        wEE=betaE*omegaEE*(1+alpha*hierNorm)
        wIE=betaI*omegaIE*(1+alpha*hierNorm)
        wEI=-betaE*omegaEI; 
        wII=-betaI*omegaII; 
        
        # Synaptic weights for inter-areal connections
        wEElr=(fln.T*(betaE*muEE*(1+alpha*hierNorm))).T
        wIElr=(fln.T*(betaI*muIE*(1+alpha*hierNorm))).T
        
        wEEaux=np.diag(-1+wEE)+wEElr
        wEIaux=wEI*np.eye(nNodes)
        
        wIEaux=np.diag(wIE)+wIElr
        wIIaux=(-1+wII)*np.eye(nNodes)
        
        # temp matrices to create matrix A
        A1 = np.concatenate((wEEaux/tauE, wEIaux/tauE),axis=1)
        A2 = np.concatenate((wIEaux/tauI, wIIaux/tauI),axis=1)
        A  = np.concatenate((A1, A2))
        
        B=A
        B[0:nNodes,:]=B[0:nNodes,:]*tauE
        B[nNodes:,:]=B[nNodes:,:]*tauI  
        
        curr=np.dot(-B,desiredSs)    
        
        # steady state
        sstate,_,_,_=np.linalg.lstsq(-B,curr,rcond=None)
        initCond=np.squeeze(sstate)
        
        popRate,_=rateModel(initCond,wEE,wEI,wIE,wII,wEElr,wIElr,nSteps,dt,noise,areaidx-1,areaInput,curr,tauE,tauI)     
    
        if u==0:        
            if muEE==34 or muEE==36:
                Area1data[idxMu,:]=popRate[0,:]
                Area2data[idxMu,:]=popRate[nNodes-1,:]    
                idxMu=idxMu+1
        
        Area2peak[u,muIdx]=np.max(popRate[nNodes-1,toplotstart:toplotend]-popRate[nNodes-1,toplotstart])
    
################################# PLOT ########################################

  
for i in range(2):        
    rateV1 = np.maximum(1e-2,Area1data[i,toplotstart:toplotend]-Area1data[i,toplotstart])
    rate24 = np.maximum(1e-2,Area2data[i,toplotstart:toplotend]-Area2data[i,toplotstart])

    ax=plt.subplot(2,2,i+1) 
    ax.semilogy(np.arange(-0.25,(len(rateV1)*dt)-0.25,dt),rateV1, 'dodgerblue') 
    ax.semilogy(np.arange(-0.25,(len(rateV1)*dt)-0.25,dt),rate24,'forestgreen')
    ax.set_ylim([1e-2,1e2+100])
    ax.set_xlim([-0.25,2.25])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    
    ax.set_ylabel('Change in firing rate  (Hz)', fontsize='large')
    ax.set_xlabel('Time (s)', fontsize='large')    
   
     
Area2peak[Area2peak>500]=500  
ax=plt.subplot(2,2,3) 
ax.semilogy(muEErange[4],np.squeeze(Area2peak[0,4]), 'cornflowerblue', marker="o",  markersize=12, markerfacecolor='w')
ax.semilogy(muEErange,np.squeeze(Area2peak[0,:]), 'cornflowerblue', marker=".",  markersize=10)
ax.semilogy(muEErange[3],np.squeeze(Area2peak[0,3]), 'cornflowerblue', marker="x",  markersize=15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([1e-6,1e3])
ax.set_ylabel('Maximum rate in 24c (Hz)', fontsize='large')
ax.set_xlabel('Global E to E coupling', fontsize='large')  

ax=plt.subplot(2,2,4)     
ax.set_title('Without feedback')
ax.semilogy(muEErange,np.squeeze(Area2peak[1,:]), 'cornflowerblue',marker=".",  markersize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([1e-6,1e3])
ax.set_ylabel('Maximum rate in 24c (Hz)', fontsize='large')
ax.set_xlabel('Global E to E coupling', fontsize='large')  

plt.subplots_adjust(hspace=0.5,wspace=0.3)


# Maximise the plotting window
plot_backend = matplotlib.get_backend()
mng = plt.get_current_fig_manager()
if plot_backend == 'TkAgg':
    mng.resize(*mng.window.maxsize())
elif plot_backend == 'wxAgg':
    mng.frame.Maximize(True)
elif plot_backend == 'Qt4Agg':
    mng.window.showMaximized()


plt.show() 
        