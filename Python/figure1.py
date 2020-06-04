import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import os.path
import matplotlib

#============================= Rate  Model ===================================#
def firingRate(A,trange,tstep):
    # A       - connectivity matrix
    # trange  - time range for simulation
    # tstep   - time step 

    init=np.array([1,0])                        # initial conditions
    state=np.zeros((2,len(trange)))
    state[:,0]=init
    tmpZeros=np.zeros((2,1))
    for t in range(len(trange)-1):
        
        state[:,t+1]=np.maximum(state[:,t] + tstep*A@state[:,t] , tmpZeros.T)                   
    
    return state

#=================== Figure 1B in the paper ==================================#

### Parameters
    
tstep = 0.0001;                                  # time step 
tmaxFig1B = 0.6                                  # time length
trangeFig1B=np.arange(0,tmaxFig1B,tstep)         # time points 
tE=0.02                                          # time constant 

wEEweak=4.45                                # synpatic weight E to E weak LBA    
wEIweak=4.7                                 # synpatic weight I to E weak LBA
wEEstrong=6                                 # synaptic weight E to E strong LBA
wEIstrong=6.7                               # synaptic weight I to E strong LBA
wIE=4 + 2.0/7.0                             # synpatic weight E to I
wII=wIE*1.1                                 # synpatic weight I to I 

# connectivity matrices
Aweak=(1/tE)*np.array([[wEEweak-1, -wEIweak],[wIE, -wII-1]])
Astrong=(1/tE)*np.array([[wEEstrong-1, -wEIstrong],[wIE, -wII-1]])

# Rate Model 
stateWeak=firingRate(Aweak,trangeFig1B,tstep)
stateStrong=firingRate(Astrong,trangeFig1B,tstep)


#=================== Figure 1C in the paper ==================================#

# Paramters

tmaxFig1C = 2                                    # time length
trange=np.arange(0,tmaxFig1C,tstep)              # time points
wEErange=np.arange(4,6.5,.025)              # range of synaptic weights E to E
wEIrange=np.arange(4.5,7,.025)              # range of synaptic weights I to E

# Arrays to store data
maxamplocal = np.zeros((len(wEIrange),len(wEErange)));  # max amplitude of E population rate  
eigA = np.zeros((len(wEIrange),len(wEErange)))          # eigenvalues

# Run simulation for range of synaptic weights

if not os.path.isfile('maxamplocal.npy'):
    for i in range(len(wEErange)):
        print(i)
        
        for j in range(len(wEIrange)):
            
            wEE=wEErange[i]
            wEI=wEIrange[j]
            
            A=np.array([[wEE-1, -wEI],[wIE, -wII-1]])
            A=(1/tE)*A
            
            state=firingRate(A,trange,tstep)
    
            wcond,x =la.eig(A);
            eigA[j,i] = np.max(np.diag(np.real(np.diag(wcond))));
    
            if eigA[j,i]<0:
                maxamplocal[j,i] = np.max(state[0,:])
    
    # set value -5 to unstable region (it helps to plot the colormap)            
    for i in range(len(wEErange)):
        for j in range(len(wEIrange)):
            if maxamplocal[j,i] ==0:  
                maxamplocal[j,i]=-5

    # Save file
    np.save('maxamplocal.npy', maxamplocal)
else:
    maxamplocal=np.load('maxamplocal.npy')

    
#=================================  Panel ===================================#

ax=plt.subplot(223)
ax.plot(trangeFig1B,stateWeak[0,:], 'g')
ax.plot(trangeFig1B,stateStrong[0,:],'m')
ax.set_ylabel('Excitatory rate (Hz)', fontsize='x-large')
ax.set_xlabel('Time (s)', fontsize='x-large')
ax.set_ylim([0,6])
ax.set_xlim([0,tmaxFig1B])
ax.set_yticks([0,2,4,6])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(['weak LBA','strong LBA'],prop={'size': 15},frameon=False)

plt.subplot(122)
X, Y = np.meshgrid(wEErange, wEIrange)
levels = np.linspace(-5,8,20)
plt.contourf(X,Y,maxamplocal,levels=levels,cmap='hot')       
plt.contourf(X,Y,maxamplocal,levels=np.linspace(-5,0),colors='silver')     
plt.plot(wEEstrong,wEIstrong,'mx')
plt.plot(wEEweak,wEIweak,'gx')
plt.ylabel('Local I to E coupling', fontsize=15)
plt.xlabel('Local E to E coupling', fontsize=15)


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