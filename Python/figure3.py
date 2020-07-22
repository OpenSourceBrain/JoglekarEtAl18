import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib
import os

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

# Path for .mat files
path=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Matlab/'

# anatomical
data=loadmat(path+'subgraphData.mat')
nNodes=data['nNodes'][0,0]
hier=np.squeeze(data['hierVals'])
hierNorm=hier/max(hier)
fln=data['flnMat']

areaList=[]
for i in range(nNodes):
    areaList.append(data['areaList'][i][0][0])

# Input 
tStart=2      
tEnd=2.25
ampl=21.8*1.9
nSteps=5e4
dt=2e-4
tSteps=np.arange(dt,dt*nSteps,dt);
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
noise=0
# Values or rate at steady state
desiredSs=np.concatenate((10*np.ones((nNodes,1)), 35*np.ones((nNodes,1)))) 

tauE=2e-2
tauI=1e-2

sStrong=1000
sWeak=100

areaidx=1


toplotstart=round(1.75/dt);
toplotend=round(5/dt);

# Store data for figure B
# (strong or weak) x nNodes x nSteps 
areaData=np.zeros((2,nNodes,int(nSteps)))
areaPeak=np.zeros((2,nNodes))

### Figure B, C and D
idxMu=0
for u in range(2):
    
    if u==0: # weak GBA
        omegaEI = 19.7 
        muEE = 33.7 
        ampl = 22.05*1.9 
        areaInput[round(tStart/dt) : round(tEnd/dt) ]=ampl
    else: # strong GBA
        omegaEI = 25.2 
        muEE = 51.5 
        ampl = 11.54*1.9 
        areaInput[round(tStart/dt) : round(tEnd/dt) ]=ampl
        
        
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
    
    areaData[u,:,:],_=rateModel(initCond,wEE,wEI,wIE,wII,wEElr,wIElr,nSteps,dt,noise,areaidx-1,areaInput,curr,tauE,tauI)     

    for i in range(nNodes):
        areaPeak[u,i]=np.max(areaData[u,i,toplotstart:toplotend]-areaData[u,i,toplotstart])


# Figure E
ampl=22.05*1.9
areaInput=np.zeros((int(nSteps)))
areaInput[round(tStart/dt) : round(tEnd/dt) ]=ampl
muEErange = np.arange(20,52,2)
omegaEI = 19.7
Area2peak=np.zeros((2,len(muEErange)))
for u in range(2):
    for muIdx in range(len(muEErange)): #len(muEErange)

        
        muEE=muEErange[muIdx]
        
        if u==1: 
            omegaEI = 19.7 + (muEE-33.7)*55/178;  
        
                 
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
    
        Area2peak[u,muIdx]=np.max(popRate[nNodes-1,toplotstart:toplotend]-popRate[nNodes-1,toplotstart])
    




##################### Plot ####################################################

# Figure B
listAreas=[0,2,5,7,8,12,16,18,28]
for i,j in enumerate(listAreas):
    ax=plt.subplot(18,2,2*i+1);
    ax.plot(tSteps[toplotstart:toplotend]-100,areaData[0,j,toplotstart:toplotend]-areaData[0,j,toplotstart],'green')
    ax.set_xlim([-98.25,-96])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)                   
    ax.spines['bottom'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    if j==0:
        ax.set_ylim([0,140])
        ax.set_yticks([round(sWeak*areaPeak[0,j],1)/sWeak])
        ax.set_title('Weak GBA')
    else:
        ax.set_ylim([0,1.2*areaPeak[0,j]])
        ax.set_yticks([round(sWeak*areaPeak[0,j],1)/sWeak])
        
    if i==4:
         ax.set_ylabel('Change in firing rate  (Hz)', fontsize='large')

# Figure C
listAreas=[0,2,5,7,8,12,16,18,28]
for i,j in enumerate(listAreas):
    ax=plt.subplot(18,2,2*i+2);
    ax.plot(tSteps[toplotstart:toplotend]-100,areaData[0,j,toplotstart:toplotend]-areaData[0,j,toplotstart],'green')
    ax.plot(tSteps[toplotstart:toplotend]-100,areaData[1,j,toplotstart:toplotend]-areaData[1,j,toplotstart], 'purple')
    ax.set_xlim([-98.25,-96])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)                   
    ax.spines['bottom'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    if j==0:
        ax.set_ylim([0,140])
        ax.set_yticks([round(sStrong*areaPeak[1,j])/sStrong])
        ax.set_title('Strong GBA')
    else:
        ax.set_ylim([0,1.2*areaPeak[1,j]])
        ax.set_yticks([round(sStrong*areaPeak[1,j])/sStrong])
    
    ax.annotate(areaList[j], xy=(-96, 1.2*areaPeak[1,j]/3))    
        
# Figure D
ax=plt.subplot(3,2,5)     
ax.semilogy(np.arange(0,nNodes),100*areaPeak[0,:]/areaPeak[0,0],'green',marker=".",  markersize=5)
ax.semilogy(np.arange(0,nNodes),100*areaPeak[1,:]/areaPeak[1,0],'purple',marker=".",  markersize=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)                    
ax.set_ylim([1e-4,1e2])
ax.set_xlim([0,nNodes])
ax.legend(['weak GBA','strong GBA'],prop={'size': 10},loc='upper right',bbox_to_anchor=(1.0, 1.2),frameon=False)
ax.set_xticks(np.arange(0,nNodes))
ax.set_xticklabels(areaList,rotation='vertical', fontsize=10)

# Figure E
Area2peak[Area2peak>500]=500  
ax=plt.subplot(2,3,6) 
ax.semilogy(muEErange,np.squeeze(Area2peak[0,:]),'cornflowerblue',marker=".",  markersize=5)
ax.semilogy(muEErange,np.squeeze(Area2peak[1,:]), 'black',marker=".",  markersize=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([1e-8,1e4])
ax.set_xlim([20,50])
ax.set_ylabel('Maximum rate in 24c (Hz)', fontsize='large')
ax.set_xlabel('Global E to E coupling', fontsize='large')  


plt.show()
