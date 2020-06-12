# -*- coding: utf-8 -*-
"""
 to create figures for spiking network models in joglekar et al neuron 2018
@author: maddy -- generate figures for paper. 
"""
from __future__ import division
from brian2 import *

import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['font.family'] = "sans-serif"

import scipy.io
import numpy as np


def cm2inch(value):
    return value/2.54
    
rate = -2 #counting... asyn to synchronous case. 

arealist = ['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c']

#spiking model - synchronous regime 
if rate == 0:
    badprop=load('/Users/maddy/Dropbox/rishicode_maddy070115/allspikeallspikebadnew2.npy')
    goodprop=load('/Users/maddy/Dropbox/rishicode_maddy070115/allspikegood.npy')
    duration = 500*ms

#spiking model - asynchronous regime 
if rate == 1:  
    badprop = load('/Users/maddy/Dropbox/rishicode_maddy070115/allspikecur15_150msbadrate.npy') 
    goodprop = load('/Users/maddy/Dropbox/rishicode_maddy070115/allspikecur6pt3_150msgoodrate.npy') 
    duration = 650*ms

#no of areas.. 
arealen = 29    
binsize, stepsize = 10*ms, 1*ms  
poster, fsize  = 0, 7

if rate > -1: 
        #####for peak firing rates plot. ##################
        maxratebad, maxrategood = np.empty([arealen,1]), np.empty([arealen,1])
        
        #sort net spikes
        netspikebad, netspikegood = len(badprop), len(goodprop)
        badpropsorted = badprop[badprop[:,1].argsort(),]
        goodpropsorted = goodprop[goodprop[:,1].argsort(),]
        
        netbinno = int( 1+(duration/ms)-(binsize/ms))
        popratebad = np.empty([arealen,netbinno ])
        poprategood = np.empty([arealen,netbinno ])
        
        u = 0 #for areas. 
        countbad, countgood = 0, 0#for each spike. 
        
        monareaktimeallbad, monareaktimeallgood = [], []
                
        while u<arealen:
            monareaktimebad, monareaktimegood = [], []
            
            while((countbad < netspikebad) and (badpropsorted[countbad,1]<1600*(u+1)) ):
              monareaktimebad.append(badpropsorted[countbad,0])#append spike times for each area.
              countbad = countbad + 1
              
            while((countgood < netspikegood) and (goodpropsorted[countgood,1]<1600*(u+1)) ):
              monareaktimegood.append(goodpropsorted[countgood,0])#append spike times for each area.
              countgood = countgood + 1
              
            valsbad, valsgood = [], []
            valsbad = numpy.histogram(monareaktimebad, bins=int(duration/stepsize))
            valsgood = numpy.histogram(monareaktimegood, bins=int(duration/stepsize))
            
            valszerobad, valszerogood = valsbad[0], valsgood[0]

            astep = binsize/(1*ms)
            valsnewbad, valsnewgood = np.zeros(netbinno), np.zeros(netbinno)
            
            acount = 0
            while acount < netbinno:        
                valsnewbad[acount] = sum(valszerobad[int(acount):int(acount+astep)])
                valsnewgood[acount] = sum(valszerogood[int(acount):int(acount+astep)])
                acount=acount+1
        
            valsratebad = valsnewbad*((1000*ms/binsize) /(1600) ) # divide by no of neurons per E pop. 
            valsrategood = valsnewgood*((1000*ms/binsize) /(1600) )    
            popratebad[u,:], poprategood[u,:] = valsratebad, valsrategood    
            #compute population firing rates. 
            
            maxratebad[u,0] = max(valsratebad[int(len(valsratebad)/3):])
            maxrategood[u,0] = max(valsrategood[int(len(valsrategood)/3):])
                
            monareaktimeallbad.append(monareaktimebad) #append for each area. 
            monareaktimeallgood.append(monareaktimegood)
            u = u+1
            
        print maxratebad[0],maxrategood[0]
        
        #"""
        badpropfig = plt.figure(figsize=(cm2inch(5.5), cm2inch(8.5)))
        ax = badpropfig.add_subplot(1,1,1)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        every, msz = 5, .05 
        
        #weak gba figure
        plt.plot(badprop[::every,0], 1.0*badprop[::every,1]/(1600), '.',markersize=msz, color = (.0,0.447,0.741))
        plt.plot([0, duration/ms], np.arange(arealen+1).repeat(2).reshape(-1, 2).T, 'k-')
        plt.ylabel('Area',fontsize=fsize)
        plt.yticks(np.arange(arealen))
        plt.xlabel('Time (ms)',fontsize=fsize)
        ylim(0,arealen)
        h=yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5],['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'],fontsize = fsize)
        cb = (241/256., 163/256., 64/256.)
        h[1][0].set_color(cb)
        h[1][1].set_color(cb)
        h[1][2].set_color(cb)
        h[1][8].set_color(cb)
        h[1][18].set_color(cb)

        if rate == 0:
            xlim(300,355)
            xticks([300, 350],[0,50],fontsize = fsize)
            plt.show()
        #badpropfig.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 5A.pdf',bbox_inches='tight', pad_inches=0.05)
        
        #strong gba figure
        goodpropfig = plt.figure(figsize=(cm2inch(5.5), cm2inch(8.5)))
        ax = goodpropfig.add_subplot(1,1,1)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')        
        plt.plot(goodprop[::every,0], 1.0*goodprop[::every,1]/(1600), '.',markersize=msz, color = (.0,0.447,0.741))
        plt.plot([0, duration/ms], np.arange(arealen+1).repeat(2).reshape(-1, 2).T, 'k-')
        plt.xlabel('Time (ms)',fontsize=fsize)
        ylim(0,arealen)
        yticks([],[],fontsize = fsize)

        if rate == 0:
            xlim(300,355)
            xticks([300, 350],[0,50],fontsize = fsize)
            plt.show()
            #goodpropfig.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/new main fig 5 goodprop.pdf',bbox_inches='tight', pad_inches=0.05)

        #peak firing rates across areas.. weak and strong gba  
        lw = 1
        atte = plt.figure(figsize=(cm2inch(7.5), cm2inch(4.)))
        ax = atte.add_subplot(1,1,1)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_yscale('log')
        plt.plot(range(arealen), maxratebad, '.-',color = (0.4660, 0.6740, 0.1880),linewidth = lw,label='weak GBA')#'Color',[0.4660, 0.6740, 0.1880]);,'Color',[0.4940,0.1840,0.5560]
        plt.plot(range(arealen), maxrategood,'.-', color = (0.4940,0.1840,0.5560),linewidth = lw,label='strong GBA')
        plt.xlabel('Area',fontsize=7)
        plt.xticks(np.arange(arealen))
        xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5]-.5*ones(arealen),['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'],rotation=90, fontsize = 7)
        yticks(fontsize = 7)
        ylim(1,250)
        xlim(0,29)
        plt.legend(fontsize = 7,frameon=False,bbox_to_anchor=(.62, 1.05))
        plt.ylabel('Maximum firing rate (Hz)', fontsize=7)
        plt.minorticks_off()
        plt.show()
        if rate == 1: #for caret figures. 
        
        #   atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 4E.pdf',bbox_inches='tight', pad_inches=0.05)
           bad = maxratebad/max(maxratebad)
           good = maxrategood/max(maxrategood)
           brate = -log10(bad)/max(-log10(bad))
           grate = -log10(good)/max(-log10(good))
#           scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/brate.mat', mdict={'brate': brate})
#           scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/grate.mat', mdict={'grate': grate})
           
        if rate == 0:    
        #   atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 5D.pdf',bbox_inches='tight', pad_inches=0.05)
        #atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 5 attevsarea.pdf',bbox_inches='tight', pad_inches=0.05)
           bad = maxratebad/max(maxratebad)
           good = maxrategood/max(maxrategood)
           bspike = -log10(bad)/max(-log10(bad))
           gspike = -log10(good)/max(-log10(good))
#           scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/bspike.mat', mdict={'bspike': bspike})
#           scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/gspike.mat', mdict={'gspike': gspike})
        
        
        #figures for 
        #onset time
        if rate == 0:
            #mean rates. 
            meanrategood = np.empty([arealen,1])
            spantimegood, spanstarttimegood = np.zeros([arealen,1]), np.zeros([arealen,1])
            
            netspikegood = len(goodprop)
            goodpropsorted = goodprop[goodprop[:,1].argsort(),]
            
            netbinno = int( 1+(duration/ms)-(binsize/ms))
                
            poprategood = np.empty([arealen,netbinno ])
            
            u = 0 #for areas. 
            countgood = 0#for each spike. 
            
            monareaktimeallgood = []
            
            while u<arealen:
                monareaktimegood =  []
                
                while((countgood < netspikegood) and (goodpropsorted[countgood,1]<1600*(u+1)) ):
                  monareaktimegood.append(goodpropsorted[countgood,0])#append spike times. for each area.
                  countgood = countgood + 1
                  
                valsgood =  []
                valsgood = numpy.histogram(monareaktimegood, bins=int(duration/stepsize))
                
                valszerogood =  valsgood[0]
            
                astep = int(binsize/stepsize)
                valsnewgood = np.zeros(netbinno)
                
                acount = 0
                while acount < netbinno:        
                    valsnewgood[acount] = sum(valszerogood[acount:acount+astep])
                    acount=acount+1
            
                valsrategood = valsnewgood*((1000*ms/binsize) /(1600) )    
                
                poprategood[u,:] =  valsrategood    
               
                #compute max and mean rates    
                maxrategood[u,0] = max(valsrategood[int(len(valsrategood)/3):])
                #meanrategood[u,0] = mean(valsrategood[int(len(valsrategood)/3):int(3*len(valsrategood)/5)])
                meanrategood[u,0] = mean(valsrategood[int(len(valsrategood)/3):int(.55*len(valsrategood))])
            
                #find index for max rate
                arrind, arrindgood = [], []
                arrindgood = np.where(poprategood[u,:]>1.5*meanrategood[u,0])
                
                maxrateindexgood = np.where(poprategood[u,:]==maxrategood[u,0])
                
                #check if highest value is at many points. 
                checkmultgood =  0
                if len(maxrateindexgood[0]) > 0:
                    for inlen in range(len(maxrateindexgood[0]) -1 ):
                        if (maxrateindexgood[0][inlen+1] - maxrateindexgood[0][inlen] > 1):
                            checkmultgood = 1
                            print inlen, "inlen"
                if checkmultgood == 0:
                    maxrateindexgoodnet = maxrateindexgood[0][0]
                if maxrateindexgoodnet < 300:
                    checkmultgood = 1    
                    print "not prop", u
                if checkmultgood == 0 : #len(maxrateindexgood[0])==1:
                    locofmaxindgood = np.where(arrindgood[0] == maxrateindexgoodnet )[0][0] 
                    #find closest indices on both sides. 
                    locnewplus, locnewminus, togp, togm = locofmaxindgood, locofmaxindgood, 0, 0
                    while (togp==0 and locnewplus<len(arrindgood[0]) -1 ):
                      locnewplus = locnewplus + 1
                      if ( arrindgood[0][locnewplus] - arrindgood[0][locofmaxindgood]> locnewplus - locofmaxindgood):
                         togp = 1
                    
                    while (togm==0 and locnewminus > 0):
                      locnewminus = locnewminus - 1
                      if ( arrindgood[0][locofmaxindgood] - arrindgood[0][locnewminus]> locofmaxindgood - locnewminus):
                         togm = 1     
                    
                    locnewplusnet = locnewplus - 1
                    locnewminusnet = locnewminus + 1
                    spantimegood[u,0] = locofmaxindgood
                    spanstarttimegood[u,0] = arrindgood[0][locnewminusnet]
            
                    print u, spanstarttimegood[u,0]
            
            #now find set of indices around max rate index lying in arrind.     
                monareaktimeallgood.append(monareaktimegood)
                u = u+1
            
            #print spantimegood
            listnotprop = []
            
            list1spike = [15,22]
            list1spike=[22]#this area doesnt show activity
            
            for u in range(29): 
             if u in list1spike:   
                   spanstarttimegood[u] = .5*(spanstarttimegood[u-1]+spanstarttimegood[u+1] )       
                   listnotprop.append(u)
             
            
            flnMatp = scipy.io.loadmat('/Users/maddy/anaconda2/efelenMatpython.mat')
            flnMat=flnMatp['flnMatpython'][:][:] #fln values..Cij is strength from j to i 
            
            flnMat = flnMat*(1.*(flnMat > .02))
            
            distMatp = scipy.io.loadmat('/Users/maddy/anaconda2/subgraphWiring29.mat')
            distMat=distMatp['wiring'][:][:] #dist values..
            delayMat = distMat/3.5
            posfln = (flnMat.T > 0)*1.
            
            delaywithfln = delayMat*posfln
            
            import scipy.sparse.csgraph
            output, predecessors = scipy.sparse.csgraph.shortest_path(delaywithfln,return_predecessors=True)
            
            shortestlen = np.zeros([29])
            for a in range(1,29):
                k=1
                achange = a
                while predecessors[0,achange] != 0:
                    achange = predecessors[0,achange]
                    k = k + 1
                shortestlen[a] = k
            #print short
            lw=1
            
            atte = plt.figure( figsize=(cm2inch(7.8), cm2inch(2.)) )#12.5 6.5
            ax = atte.add_subplot(1,1,1)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            plt.plot(range(29),shortestlen + output[0,:],'.-',linewidth = lw,color = (241/256., 163/256., 64/256.),label='predicted')# (FLN threshold 0.02)')
            plt.plot(range(29),spanstarttimegood-spanstarttimegood[0],'.-',linewidth = lw,color = (0.4940,0.1840,0.5560),label='observed')
            
            for uu in list1spike: #listnotprop:
              #print uu  
              plt.plot(range(uu, uu+1),spanstarttimegood[uu:1+uu]-spanstarttimegood[0],'.r')
            plt.xlabel('Area',fontsize=7)
            plt.xticks(np.arange(arealen))
            xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5]-.5*ones(arealen),['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'],rotation=90, fontsize = 7)
            yticks([0, 30, 60],[0, 30, 60],fontsize = 7)
            plt.ylabel('Onset time (ms)', fontsize=7)
            xlim(0,29)
            
            plt.minorticks_off()
            plt.legend(fontsize = 7,frameon=False,bbox_to_anchor=(.38, 1.35))#bbox_to_anchor=(.3, 1.07)
            plt.show()
            atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 5 onset1.pdf',bbox_inches='tight', pad_inches=0.05)

            atte = plt.figure( figsize=(cm2inch(3), cm2inch(3)) ) #1.75 1.75
            ax = atte.add_subplot(1,1,1)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            plt.plot(spanstarttimegood-spanstarttimegood[0],spanstarttimegood-spanstarttimegood[0],'-',linewidth = lw,color = 'y')
            #plt.plot(spanstarttimegood-spanstarttimegood[0],1.2*output[0,:],'.',color = (241/256., 163/256., 64/256.))# (FLN threshold 0.02)')
            plt.plot(spanstarttimegood-spanstarttimegood[0],shortestlen+output[0,:],'.',linewidth = lw,color = (241/256., 163/256., 64/256.))# (FLN threshold 0.02)')
            
            for uu in list1spike: #listnotprop:
                  #plt.plot(spanstarttimegood[uu:1+uu]-spanstarttimegood[0],output[0,u],'.r')
                  plt.plot(spanstarttimegood[uu:1+uu]-spanstarttimegood[0],shortestlen[u]+output[0,u],'.r')
            
            plt.xlabel('Observed onset (ms)',fontsize=7)
            plt.ylabel('Predicted onset (ms)',fontsize=7)
            #plt.title('Onset time (ms)',fontsize = 7)
            xticks([0, 30, 60],[0, 30, 60],fontsize = 7)
            yticks([0, 30, 60],[0, 30, 60],fontsize = 7)

            plt.minorticks_off()
            plt.show()
                    
            sumofsq = 0
            for a in range(29):
                if a not in list1spike:
                    sumofsq = sumofsq + (output[0,a] - (spanstarttimegood[a:1+a]-spanstarttimegood[0]) )**2
            print np.sqrt(sumofsq)
            #atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 5 onset2.pdf',bbox_inches='tight', pad_inches=0.05)

            #new added strong fln path
            #"""
            flnMatp = scipy.io.loadmat('/Users/maddy/anaconda2/efelenMatpython.mat')
            flnMat=flnMatp['flnMatpython'][:][:] #fln values..Cij is strength from j to i 
            
            posfln = (flnMat.T > 0)*1.
            delaywithfln = delayMat*posfln

            flnMatgood = flnMat*(1.*(flnMat > .02))
            posflngood = (flnMatgood.T > 0)*1.
            delaywithflngood = delayMat*posflngood
            
            flnMatnew = np.zeros([29,29])
            for i in range(29):
                for j in range(29):

                    if flnMat[i,j] > 0:    
                        flnMatnew[i,j] = 1./flnMat[i,j]
            
            f = flnMatnew.T #look only at stronget paths. forget distance. 
            f = f.copy(order='C') #set flags -- check this 
            
            #various methods of predicted vs actual onset time comparison
            output, predecessors = scipy.sparse.csgraph.shortest_path(f,return_predecessors=True)
            outputgood, predecessorsgood = scipy.sparse.csgraph.shortest_path(delaywithflngood,return_predecessors=True)
            outputallfln, predecessorsallfln = scipy.sparse.csgraph.shortest_path(delaywithfln,return_predecessors=True)

            shortestlenallfln = np.zeros([29])
            for a in range(1,29):
                k=1
                achange = a
                while predecessorsallfln[0,achange] != 0:
                    achange = predecessorsallfln[0,achange]
                    k = k + 1
                shortestlenallfln[a] = k
                              
            shortestlengood = np.zeros([29])
            for a in range(1,29):
                k=1
                achange = a
                while predecessorsgood[0,achange] != 0:
                    achange = predecessorsgood[0,achange]
                    k = k + 1
                shortestlengood[a] = k                     
            
            shortestlen = np.zeros([29])
            pathlen = np.zeros([29])
            for a in range(1,29):
                k=1
                achange = a
                while predecessors[0,achange] != 0:
                    pathlen[a] = pathlen[a] + delaywithfln[predecessors[0,achange], achange]
                    achange = predecessors[0,achange]
                    k = k + 1
                pathlen[a] = pathlen[a] + delaywithfln[predecessors[0,achange], achange]    
                shortestlen[a] = k

            atte = plt.figure( figsize=(cm2inch(9.5), cm2inch(7.5)) )#12.5 6.5
            ax = atte.add_subplot(1,1,1)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

            #compare methods
            plt.plot(range(29),shortestlenallfln + outputallfln[0,:],'.-',linewidth=1,color = 'g',label='predicted: shortest path')# (FLN threshold 0.02)')
            plt.plot(range(29),shortestlengood + outputgood[0,:],'.-',linewidth=1,color = (241/256., 163/256., 64/256.),label='predicted: shortest path (threshold FLNs)')# (FLN threshold 0.02)')
            plt.plot(range(29),shortestlen + pathlen,'.-',linewidth=1,color = 'k',label='predicted: strongest FLNs')# (FLN threshold 0.02)')
            plt.plot(range(29),spanstarttimegood-spanstarttimegood[0],'.-',linewidth=1,color = (0.4940,0.1840,0.5560),label='observed')
            
            for uu in list1spike: #listnotprop:
              #print uu  
              plt.plot(range(uu, uu+1),spanstarttimegood[uu:1+uu]-spanstarttimegood[0],'.r')
            plt.xlabel('Area',fontsize=7)
            plt.xticks(np.arange(arealen))
            xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5]-.5*ones(arealen),['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'],rotation=90, fontsize = 7)
            yticks([0, 30, 60],[0, 30, 60],fontsize = 7)
            plt.ylabel('Onset time (ms)', fontsize=7)
            xlim(0,29)
            
            plt.minorticks_off()
            plt.legend(fontsize = 7,frameon=False, bbox_to_anchor=(.85, 1.16))
            plt.show()
            #atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/suppfig6.pdf',bbox_inches='tight', pad_inches=0.05)
            #"""
        
if rate == -1:
        """
        #computing xi. going from asynchronous to synchronous regime
        #below are values used while running simulation 
        #note for asynchronous regime 
        #para['muIEsp'],para['omegaIIsp'],para['omegaEEsp'],para['omegaIEsp'] =.19/4*mV, .075*mV,.01*mV, .075*mV
        #below asy: asynchronous and sync: synchronous
        
        para['muIEspasy'], para['omegaIIspasy'], para['omegaEEspasy'], para['omegaIEspasy'] = para['muIEsp'], para['omegaIIsp'], para['omegaEEsp'], para['omegaIEsp']
        para['omegaEIspasy'], para['muEEspasy'] = .0375*mV, .0375*mV
        para['omegaEIspsync'], para['muEEspsync'] = .56*mV, .16*mV
        para['muEasy'], para['muIasy'] =  14.2*mV, 14.7*mV
        para['muEsync'], para['muIsync'] =  15.4*mV, 14.*mV
        currvalasy,currdurasy = 15., 1500
        currvalsync,currdursync = 10., 1000
        steps = 1.5  counting = .25 
        scal = 1. + steps*counting
        para['muIEsp'],para['omegaIIsp'],para['omegaEEsp'],para['omegaIEsp'] = scal*para['muIEspasy'], scal*para['omegaIIspasy'], scal*para['omegaEEspasy'], scal*para['omegaIEspasy']
        para['omegaEIsp'] = para['omegaEIspasy'] + (scal-1)*(1/3)*(para['omegaEIspsync'] - para['omegaEIspasy'])
        para['muEEsp'] = para['muEEspasy'] + (scal-1)*(1/3)*(para['muEEspsync'] - para['muEEspasy'])
        currval = currvalasy + (scal-1)*(1/3)*(currvalsync - currvalasy) 
        currdur = currdurasy + (scal-1)*(1/3)*(currdursync - currdurasy)
        para['muE'] = para['muEasy'] + (scal-1)*(1/3)*(para['muEsync'] - para['muEasy'])
        para['muI'] = para['muIasy'] + (scal-1)*(1/3)*(para['muIsync'] - para['muIasy'])

        #this gives values. 
        #rnd_seed = 1, 2, 3, 4, 5 duration 500ms 
        #counting 0    0.102570695745 0.108266970742 0.107345017765 0.101813213807 0.110493412191
        #counting .25  0.092943427808  0.10165399973 0.0989555004113 0.0930852453716 0.0989148643989
        #counting .5   0.104856480962  0.103376307306  0.109652868024 0.10129692177 0.107553409523
        #counting .75  0.139407341265 0.147645326908  0.144404545499  0.165579568531 0.136254393399
        #counting 1    0.285324986861 0.264488053586 0.315449285306 0.263747134019 .2231977553
        #counting 1.25 0.412162518196  0.380991126065  0.462325851462  0.42252185704  0.361195661674
        #counting 1.5  0.49249508386 0.482918744538  0.477361173603 0.532023636421 0.518135136951
        #counting 1.75 0.532882902449 0.560494600423 0.62464798306 0.598130702026  0.565986354493
        #counting 2    0.571971086045 0.598730542685  0.588065654116 0.685067550863  0.624081332854
        #countingsave[0,:] = [0.102570695745, 0.108266970742, 0.107345017765, 0.101813213807, 0.110493412191] and so on.
        #np.save('/Users/maddy/Dropbox/rishicode_maddy070115/countingsave',countingsave)
        """
        
        countingsave = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/countingsave.npy')
        lw=1
#        atte = plt.figure(figsize=(cm2inch(7.5), cm2inch(4.5))) 7.5 6.5
        atte = plt.figure(figsize=(cm2inch(4.7), cm2inch(2.1)))
        ax = atte.add_subplot(1,1,1)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        e = np.sqrt(var(countingsave,axis = 1))
        plt.errorbar(np.arange(0,1.1,.125), sum(countingsave, axis=1)/5, e, linestyle='-', linewidth=lw,marker='.',markersize=3,color=(0.4660, 0.6740, 0.1880),capsize=2)
            xlim(-.05,1.05)
        plt.xticks([0,1], fontsize = 7)
        xticks([0, 1],['Asynchronous\n model','Synchronous\n model'],fontsize = 7)  
        ylim(0,.67)
        yticks([.0, .3, .6],[.0, .3, .6],fontsize = 7)
        plt.ylabel('Degree of population\n synchrony', fontsize=7)
        plt.xlabel('Parameter space', fontsize=7)
        #plt.xlabel('From the Asynchronous \nto Synchronous state', fontsize=7)# \n along the parameter space
        plt.show()
#        atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/suppfig5.pdf',bbox_inches='tight', pad_inches=0.05)
        


if rate == -2:
    #for 16 active areas consciousness figure -- vary input strength
    diffbkandpeakratecurorig = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/diffratecur_clust_orig.npy')
    diffbkandpeakratecur = diffbkandpeakratecurorig[:,1:]
    #areas not showing much activity -- asyn regime
    list1rate = ['6','9','10','11','15','17','19','22','23','24','25','27','28']
    #areas not showing much activity -- sync regime
    list1spike = ['10','15','22','28']
    meanend = 285
    alp = 0.35
    listprefront, listtemp, listpar, listocc = [], [], [], []
    netprefront, nettemp, netpar, netocc = 0, 0, 0, 0
    collist = ['k','k','k','g','b','r','g','r','b','g','c','b','g','r','r','r','r','c','b','b','g','g','c','b','c','c','r','b','y'] 
    xc = np.arange(3.5,6.1,.5/3) 
    
    lw=1.
    
    atte = plt.figure(figsize=(cm2inch(12.5), cm2inch(4.5)))#6.8 2.9 8.5 7.5 ..11,7
    #    atte = plt.figure(figsize=(cm2inch(8.5), cm2inch(5.5)))#8.5 7.5 

    ax = atte.add_subplot(1,1,1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for a in np.arange(0,29):
        if str(a) not in list1rate:
                plt.plot(xc,diffbkandpeakratecur[a,:],'-',color = collist[a],alpha=alp, linewidth=lw)                    
                if collist[a] == 'r':
                  listprefront.append(diffbkandpeakratecur[a,:])
                if collist[a] == 'k':
                  listocc.append(diffbkandpeakratecur[a,:])
                if collist[a] == 'g':
                  listpar.append(diffbkandpeakratecur[a,:])
                if collist[a] == 'b':
                  listtemp.append(diffbkandpeakratecur[a,:])

    for i in range(len(listprefront)):
        netprefront = netprefront + listprefront[i]
    netprefront = netprefront/len(listprefront)
    plt.plot(xc,netprefront,'-',color = 'r',label='frontal',linewidth=lw+1)

    for i in range(len(listtemp)):
        nettemp = nettemp + listtemp[i]
    nettemp = nettemp/len(listtemp)
    plt.plot(xc,nettemp,'-',color = 'b',label='temporal',linewidth=lw+1)
    
    for i in range(len(listpar)):
        netpar = netpar + listpar[i]
    netpar = netpar/len(listpar)
    plt.plot(xc,netpar,'-',color = 'g',label='parietal',linewidth=lw+1)
    
    for i in range(len(listocc)):
        netocc = netocc + listocc[i]
    netocc = netocc/len(listocc)
    plt.plot(xc,netocc,'-',color = 'k',label='occipital',linewidth=lw+1)

    #plt.xlabel('Input strength (pA)',fontsize=7)
    xticks([3.5, 4., 4.5, 5., 5.5, 6],[70, '',90,'','', 120],fontsize = 7)

    plt.ylabel('Peak response', fontsize=7)
    yticks([0, 0.5, 1],[0, 0.5, 1],fontsize = 7)

    xticks([],[],fontsize = 7)
    ylim(-.01,1.03)  
    xlim(3.5, 6.1)
    plt.legend(fontsize = 7,frameon=False,bbox_to_anchor=(.02, .55))        
    
    plt.minorticks_off()
    plt.show()
#    atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 6A.pdf',bbox_inches='tight', pad_inches=0.05)

    #caret figures    
    badcons_4pt5 = 1-diffbkandpeakratecur[:,7].reshape(29,1)
    goodcons_5pt5 = 1-diffbkandpeakratecur[:,13].reshape(29,1)
    cons_5 = 1-diffbkandpeakratecur[:,10].reshape(29,1)
#    scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/badcons_4pt5.mat', mdict={'badcons_4pt5': badcons_4pt5})
#    scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/goodcons_5pt5.mat', mdict={'goodcons_5pt5': goodcons_5pt5})
#    scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/cons_5.mat', mdict={'cons_5': cons_5})


    #test how no of active areas depends on input current strength 
    dd = np.zeros([16,16])
    count = 0
    for a in np.arange(0,29):
        if str(a) not in list1rate:
            #print arealist[a]
            dd[count,:] = diffbkandpeakratecur[a,:]
            count+=1
    d = dd > .15 #.25#     OR .2   
    ddd = d.sum(axis=0)
    #plt.plot(ddd,'.-')
    print [i for i in ddd[1:] - ddd[:-1]]
    print [i for i in ddd]
    print [i for i in d[:,9]]
    ct = 0
    dict16ar = {}
    for a in np.arange(0,29):
        if str(a) not in list1rate:
            print ct, a, arealist[a]
            dict16ar[ct] = arealist[a]
            ct += 1
    for i in range(16):  
        print i, [dict16ar[k] for k in np.where(d[:,i] > 0)[0]]

    #occlist = [1,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3]
    #parlist = [0,0,0,0,0,0,1,1,4,4,4,4,4,4,4,4]
    #templist = [0,0,0,0,0,0,0,2,3,3,3,3,3,3,3,3]
    #prelist = [0,0,0,0,0,0,0,0,1,4,6,6,6,6,6,6]

    atte = plt.figure(figsize=(cm2inch(5.5), cm2inch(5.5)))
    ax = atte.add_subplot(1,1,1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.plot(xc,ddd,'.-', linewidth=lw)                    
    xticks([3.5, 4., 4.5, 5., 5.5, 6],[70, '',90,'','', 120],fontsize = 7)
    plt.ylabel('No. of areas activated', fontsize=7)
    plt.xlabel('Input strength (pA)',fontsize=7)
    yticks([0, 4, 8, 12, 16],[0, 4, 8, 12, 16],fontsize = 7)
    xlim(3.5, 6.1)    
    plt.minorticks_off()
    plt.show()
    #atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/testignition.pdf',bbox_inches='tight', pad_inches=0.05)
    
    plt.figure()      
    plt.plot(xc,occlist,'k.-',linewidth=lw)     
    plt.plot(xc,templist,'b.-',linewidth=lw)     
    plt.plot(xc,parlist,'g.-', linewidth=lw)     
    plt.plot(xc,prelist,'r.-',linewidth=lw)    
    #plt.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/ignitionareas.pdf',bbox_inches='tight', pad_inches=0.05)

if rate == -3:
    #consciousness figure no feedback.     
    diffratecur_clust_orig_nofb = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/diffratecur_clust_orig_nofb.npy')
    diffbkandpeakratecur = diffratecur_clust_orig_nofb[:,1:]
    #areas not propagating
    #async. regime case
    list1rate = ['6','9','10','11','15','17','19','22','23','24','25','27','28']
    list1rate_nofb = ['5','6','7','9','10','11','12','13','14','15','16','17','19','20','21','22','23','24','25','26','27','28']
    #areas not propagating sync regime case 
    list1spike = ['10','15','22','28']

    meanend = 285
    alp = 0.35
    listprefront, listtemp, listpar, listocc = [], [], [], []
    netprefront, nettemp, netpar, netocc = 0, 0, 0, 0    
    hfont = {'fontname':'Helvetica'}

    collist = ['k','k','k','g','b','r','g','r','b','g','c','b','g','r','r','r','r','c','b','b','g','g','c','b','c','c','r','b','y'] 
    xc = np.arange(3.5,6.1,.5/3) 
    atte = plt.figure(figsize=(cm2inch(12.5), cm2inch(4.8)))#8.5 7.5   8. 3.8
    ax = atte.add_subplot(1,1,1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    lw=1
    
    for a in np.arange(0,29):
        if str(a) not in list1rate:
          if str(a) not in list1rate_nofb:     
                plt.plot(xc,diffbkandpeakratecur[a,:],'-',color = collist[a],alpha = alp,linewidth=lw)                    

                if collist[a] == 'r':
                  listprefront.append(diffbkandpeakratecur[a,:])
                if collist[a] == 'k':
                  listocc.append(diffbkandpeakratecur[a,:])
                if collist[a] == 'g':
                  listpar.append(diffbkandpeakratecur[a,:])
                if collist[a] == 'b':
                  listtemp.append(diffbkandpeakratecur[a,:])   
                  
    for i in range(len(listprefront)):
        netprefront = netprefront + listprefront[i]
    if len(listprefront)>0:
        netprefront = netprefront/len(listprefront)
        plt.plot(xc,netprefront,'-',color = 'r',label='prefrontal',linewidth=lw+1)

    for i in range(len(listtemp)):
        nettemp = nettemp + listtemp[i]
    nettemp = nettemp/len(listtemp)
    plt.plot(xc,nettemp,'-',color = 'b',label='temporal',linewidth=lw+1)
    
    for i in range(len(listpar)):
        netpar = netpar + listpar[i]
    netpar = netpar/len(listpar)
    plt.plot(xc,netpar,'-',color = 'g',label='parietal',linewidth=lw+1)
    
    for i in range(len(listocc)):
        netocc = netocc + listocc[i]
    netocc = netocc/len(listocc)
    plt.plot(xc,netocc,'-',color = 'k',label='occipital',linewidth=lw+1)
    #plt.xlabel('Input strength (mA)',fontsize=7)
    plt.ylabel('Peak response', fontsize=7)
    yticks([0, 0.5, 1],[0, 0.5, 1],fontsize = 7)
    xticks([],[],fontsize = 7)
    xlim(3.5,6.1)
    ylim(-.01,1.03)   
    plt.minorticks_off()
    plt.show()    
#    atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 6B.pdf',bbox_inches='tight', pad_inches=0.05)

    nofbgoodcons_5pt5 = 1-diffbkandpeakratecur[:,13].reshape(29,1)#for caret figure
#    scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/nofbgoodcons_5pt5.mat', mdict={'nofbgoodcons_5pt5': nofbgoodcons_5pt5})


if rate == -5:
    #shuffle flns. 
    diffbkandpeakratecurflnshuf2 = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/diffratecur_clust_orig_flnshuf2.npy')
    diffbkandpeakratecur = diffbkandpeakratecurflnshuf2[:,1:]
    list1rate = ['6','9','10','11','15','17','19','22','23','24','25','27','28']    
    list1spike = ['10','15','22','28']

    #areas not propagating below     
    list1rate = ['1','3','6','22','23']
    alp = 0.25
    listprefront, listtemp, listpar, listocc = [], [], [], []
    netprefront, nettemp, netpar, netocc = 0, 0, 0, 0    
     
    meanend = 285
    lw=1
    
    collist = ['k','k','k','g','b','r','g','r','b','g','c','b','g','r','r','r','r','c','b','b','g','g','c','b','c','c','r','b','y'] 
    xc = np.arange(3.5,6.1,.5/3) 
    atte = plt.figure(figsize=(cm2inch(12.5), cm2inch(4.5)))#8.5 7.5 -- 8, 3.8
    ax = atte.add_subplot(1,1,1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for a in np.arange(0,29):
        if str(a) not in list1rate:
            plt.plot(xc,diffbkandpeakratecur[a,:],'-',color = collist[a],alpha = alp,linewidth=lw)                    

            if collist[a] == 'r':
              listprefront.append(diffbkandpeakratecur[a,:])
            if collist[a] == 'k':
              listocc.append(diffbkandpeakratecur[a,:])
            if collist[a] == 'g':
              listpar.append(diffbkandpeakratecur[a,:])
            if collist[a] == 'b':
              listtemp.append(diffbkandpeakratecur[a,:])                       
         
    for i in range(len(listprefront)):
        netprefront = netprefront + listprefront[i]
    netprefront = netprefront/len(listprefront)
    plt.plot(xc,netprefront,'-',color = 'r',label='prefrontal',linewidth=lw+1)

    for i in range(len(listtemp)):
        nettemp = nettemp + listtemp[i]
    nettemp = nettemp/len(listtemp)
    plt.plot(xc,nettemp,'-',color = 'b',label='temporal',linewidth=lw+1)
    
    for i in range(len(listpar)):
        netpar = netpar + listpar[i]
    netpar = netpar/len(listpar)
    plt.plot(xc,netpar,'-',color = 'g',label='parietal',linewidth=lw+1)
    
    for i in range(len(listocc)):
        netocc = netocc + listocc[i]
    netocc = netocc/len(listocc)
    plt.plot(xc,netocc,'-',color = 'k',label='occipital',linewidth=lw+1)

    #plt.xlabel('Input strength (mA)',fontsize=7)
    plt.xlabel('Input strength (pA)',fontsize=7)
    plt.ylabel('Peak response', fontsize=7)
    yticks([0, 0.5, 1],[0, 0.5, 1],fontsize = 7)
    xticks([3.5, 4.5, 6],[70, 90, 120],fontsize = 7)
    xticks([3.5, 4., 4.5, 5., 5.5, 6],[70, '',90,'','', 120],fontsize = 7)
    ylim(-.01,1.03)   
    xlim(3.5,6.1)

    plt.minorticks_off()
    plt.show()
#    atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/main fig 6C.pdf',bbox_inches='tight', pad_inches=0.05)


if rate == -6:
    #input to primary somatosensory cortex: area 2. 
    duration = 550*ms
    
    seedlist = [1,2]
    currbad, currgood = 6.9, 7.0
    
    for seed in seedlist: 
        if seed == 1:
            allspikebadS1 = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/allspike_S1inp_'+str(currbad)+'_seed_'+str(seed)+'.npy')
            popratebadS1 = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/poprate_curr_S1inp_'+str(currbad)+'_seed_'+str(seed)+'.npy')
        else:
            allspikegoodS1 = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/allspike_S1inp_'+str(currgood)+'_seed_'+str(seed)+'.npy')
            poprategoodS1 = np.load('/Users/maddy/Dropbox/rishicode_maddy070115/poprate_curr_S1inp_'+str(currgood)+'_seed_'+str(seed)+'.npy')
    
    maxratebadS1, maxrategoodS1 = np.zeros([29,1]), np.zeros([29,1])
    for u in range(29):
        maxratebadS1[u,0] = max(popratebadS1[u,int(len(popratebadS1[0,:])/3):])
        maxrategoodS1[u,0] = max(poprategoodS1[u,int(len(poprategoodS1[0,:])/3):])
    
    print (max(maxratebadS1), max(maxrategoodS1) )
    every, msz, lw, fsize = 1, .16, 1, 7
    
    #rasters
    plt.figure(figsize=(cm2inch(5.5), cm2inch(8.5)))
    plt.plot(allspikebadS1[::every,0], allspikebadS1[::every,1]/1600, '.',markersize=msz)
    plt.plot([0, duration/ms], np.arange(arealen+1).repeat(2).reshape(-1, 2).T, 'k-')
    plt.ylabel('Area')
    plt.yticks(np.arange(arealen))
    plt.xlabel('time (ms)')
    ylim(0,28)
    yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5],['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'],fontsize = fsize)
    xlim(300,duration/ms)
    plt.show()  
    
    plt.figure(figsize=(cm2inch(5.5), cm2inch(8.5)))
    plt.plot(allspikegoodS1[::every,0], allspikegoodS1[::every,1]/1600, '.',markersize=msz)
    plt.plot([0, duration/ms], np.arange(arealen+1).repeat(2).reshape(-1, 2).T, 'k-')
    plt.ylabel('Area')
    plt.yticks(np.arange(arealen))
    plt.xlabel('time (ms)')
    ylim(0,28)
    yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5],['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'],fontsize = fsize)
    xlim(300,duration/ms)
    plt.show() 
    
    #firing rates across areas -- weak and strong gba
    atte = plt.figure(figsize=(cm2inch(7.5), cm2inch(4.)))#7.5 4.5
    ax = atte.add_subplot(1,1,1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_yscale('log')
    plt.plot(range(arealen), maxratebadS1, '.-',linewidth = lw,color = (0.4660, 0.6740, 0.1880),label='weak GBA')#'Color',[0.4660, 0.6740, 0.1880]);,'Color',[0.4940,0.1840,0.5560]
    plt.plot(range(arealen), maxrategoodS1,'.-',linewidth = lw, color = (0.4940,0.1840,0.5560),label='strong GBA')
    plt.xlabel('Area',fontsize=7)
    plt.xticks(np.arange(arealen))
    xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5]-.5*ones(arealen),['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'],rotation=90, fontsize = 7)
    yticks(fontsize = 7)
    ylim(1,100)
    xlim(0,29)
    plt.legend(fontsize = 7,frameon=False,bbox_to_anchor=(1.05, 1.05))
    plt.ylabel('Maximum firing rate (Hz)', fontsize=7)
    plt.minorticks_off()
    plt.show()
    #atte.savefig('/Users/maddy/Dropbox/rishicode_maddy070115/signal propagation paper/paper figures v2/suppfig4.pdf',bbox_inches='tight', pad_inches=0.05)
    
    #for caret figure
    badS1 = maxratebadS1/max(maxratebadS1)
    goodS1 = maxrategoodS1/max(maxrategoodS1)
    brateS1 = -log10(badS1)/max(-log10(badS1))
    grateS1 = -log10(goodS1)/max(-log10(goodS1))
    #scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/brateS1.mat', mdict={'brateS1': brateS1})
    #scipy.io.savemat('/Users/maddy/Dropbox/rishicode_maddy070115/caret2/grateS1.mat', mdict={'grateS1': grateS1})

    
    
"""#check active areas as input is increased gradually. 
dd = np.zeros([16,16])
count = 0
for a in np.arange(0,29):
    if str(a) not in list1rate:
        #print arealist[a]
        dd[count,:] = diffbkandpeakratecur[a,:]
        count+=1
d = dd > .25     OR .2   
ddd = d.sum(axis=0)
plt.plot(ddd,'.-')
print [i for i in ddd[1:] - ddd[:-1]]
print [i for i in ddd]
print [i for i in d[:,9]]
ct = 0
dict16ar = {}
for a in np.arange(0,29):
    if str(a) not in list1rate:
        print ct, a, arealist[a]
        dict16ar[ct] = arealist[a]
        ct += 1
for i in range(16):  
    print i, [dict16ar[k] for k in np.where(d[:,i] > 0)[0]]
"""        
    
    
    