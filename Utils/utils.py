import numpy as np
from brian2 import *
import matplotlib.pyplot as plt

# Raster Plot
def rasterPlot(xValues,yValues,duration,figure,N,saveFigure,path):
    
    ticks=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,
       8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,
       16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,
       24.5,25.5,26.5,27.5,28.5]   

    areasName=['V1','V2','V4','DP','MT',
           '8m','5','8l','TEO','2','F1',
           'STPc','7A','46d','10','9/46v',
           '9/46d','F5','TEpd','PBr','7m','7B',
           'F2','STPi','PROm','F7','8B','STPr','24c']   

    plt.figure()
    plt.plot(xValues, 1.0*yValues/(4*400), '.',markersize=1)
    plt.plot([0, duration], np.arange(N+1).repeat(2).reshape(-1, 2).T, 'k-')
    plt.ylabel('Area')
    plt.yticks(np.arange(N))
    plt.xlabel('time (ms)')
    plt.ylim(0,N)
    plt.yticks(ticks[:N],areasName[:N])
    plt.xlim(0,duration)

    # Save figure
    if saveFigure== 'yes':
        plt.savefig(path+'/figures/figure'+figure+'_'+str(N)+'areas.png')     

    plt.show()

    return 0

# Plot for Maximum firing rate
def firingRatePlot(maxrategood,maxratebad,figure,N,saveFigure,path):
    
    areasName=['V1','V2','V4','DP','MT',
           '8m','5','8l','TEO','2','F1',
           'STPc','7A','46d','10','9/46v',
           '9/46d','F5','TEpd','PBr','7m','7B',
           'F2','STPi','PROm','F7','8B','STPr','24c']   
    
    
    plt.semilogy(range(N), maxratebad, '.-',color = (0.4660, 0.6740, 0.1880),linewidth = 1,label='weak GBA')
    plt.semilogy(range(N), maxrategood,'.-', color = (0.4940,0.1840,0.5560),linewidth = 1,label='strong GBA')
    plt.xticks(range(N), areasName)
    plt.xticks(rotation=90) 
    plt.xlabel('Area')
    plt.ylabel('Maximum Firing Rate (Hz)')
    plt.legend(('weak GBA', 'strong GBA'))
    
    # Save figure
    if saveFigure== 'yes':
        plt.savefig(path+'/figures/figure'+figure+'_'+str(N)+'areas.png')    
    
    plt.show()    
    
    return 0

def firingRate(N,goodprop,duration):

    binsize = 10 #[ms]
    stepsize =  1 #[ms] 
    
    # Store maximum firing rate for each area
    maxrategood = np.empty([N,1])
    
    # sort net spikes
    netspikegood= len(goodprop)
    
    goodpropsorted = goodprop[goodprop[:,1].argsort(),]
            
    netbinno = int( 1+(duration)-(binsize))
    poprategood = np.empty([N,netbinno ])
            
     
    countgood = 0#for each spike. 
            
    monareaktimeallgood = []
                    
    for u in range(N):
        monareaktimegood = []
          
        while((countgood < netspikegood) and (goodpropsorted[countgood,1]<1600*(u+1)) ):
          monareaktimegood.append(goodpropsorted[countgood,0])#append spike times for each area.
          countgood = countgood + 1
          
        valsgood = np.histogram(monareaktimegood, bins=int(duration/stepsize))
        
        valszerogood = valsgood[0]
    
        astep = binsize
        
        valsnewgood = np.zeros(netbinno)
        
        acount = 0
        while acount < netbinno:        
            valsnewgood[acount] = sum(valszerogood[int(acount):int(acount+astep)])
            acount=acount+1
    
        valsrategood = valsnewgood*((1000/binsize) /(1600) )    
        poprategood[u,:] = valsrategood    
        
        #compute population firing rates. 
        maxrategood[u,0] = max(valsrategood[int(len(valsrategood)/3):])
            
    return maxrategood,poprategood

