from neuromllite.NetworkGenerator import check_to_generate_or_run
from neuromllite import Simulation

from neuromllite import Network, Population, Projection, Cell, Synapse, InputSource, Input
from neuromllite import RandomConnectivity,RectangularRegion, RelativeLayout

import sys
import numpy as np
import shutil
import glob
import scipy.linalg as la
import os

# This function generates the overview of the network using neuromllite
def internal_connections(pops,W,syns):
    for pre in pops:
        for post in pops:

            weight = W[pops.index(post)][pops.index(pre)]
            print('Connection %s -> %s weight: %s'%(pre.id,
            post.id, weight))
            if weight!=0:

                net.projections.append(Projection(id='proj_%s_%s'%(pre.id,post.id),
                                                    presynaptic=pre.id,
                                                    postsynaptic=post.id,
                                                    synapse=syns[pre.id],
                                                    type='continuousProjection',
                                                    delay=0,
                                                    weight=weight,
                                                    random_connectivity=RandomConnectivity(probability=1)))
                


wEErange=np.arange(4,6.5,.05)              # range of synaptic weights E to E
wEIrange=np.arange(4.5,7,.1)              # range of synaptic weights I to E
wIE=4 + 2.0/7.0                             # synpatic weight E to I
wII=wIE*1.1                                 # synpatic weight I to I 
tE=0.02             


# Arrays to store data
maxamplocal = np.zeros((len(wEIrange),len(wEErange)));  # max amplitude of E population rate  
eigA = np.zeros((len(wEIrange),len(wEErange)))          # eigenvalues


for i in range(len(wEErange)):
    for j in range(len(wEIrange)):
    
         

        # Build the network
        net = Network(id='Joglekar_figure1c')
        net.notes = 'A simple rate model with E and I populations'
        
        
        r1 = RectangularRegion(id='Joglekar', x=0,y=0,z=0,width=1000,height=100,depth=1000)
        net.regions.append(r1)
       
        wEE=wEErange[i]
        wEI=wEIrange[j]  
        
        
        A=(1/tE)*np.array([[wEE-1, -wEI],[wIE, -wII-1]])
        
        
        wcond,x =la.eig(A);
        eigA[j,i] = np.max(np.diag(np.real(np.diag(wcond))));
    
        if eigA[j,i]<0:
        
            net.parameters = { 'wee':      A[0,0],
                               'wei':      A[1,0],
                               'wie':      A[0,1],
                               'wii':      A[1,1]}  
            
            
            exc_cell = Cell(id='Exc', lems_source_file='figure1c_Parameters.xml')
            inh_cell = Cell(id='Inh', lems_source_file='figure1c_Parameters.xml')
            net.cells.append(exc_cell)
            net.cells.append(inh_cell)
            
            exc_pop = Population(id='Excitatory_'+str(round(wEE,5))+'_'+str(round(wEI,5)), 
                                 size=1, 
                                 component=exc_cell.id, 
                                 properties={'color': '0.8 0 0','radius':10},
                                 relative_layout = RelativeLayout(region=r1.id,x=-20,y=0,z=0))
            
            inh_pop = Population(id='Inhibitory_'+str(round(wEE,5))+'_'+str(round(wEI,5)), 
                                 size=1, 
                                 component=inh_cell.id, 
                                 properties={'color': '0 0 0.8','radius':10},
                                 relative_layout = RelativeLayout(region=r1.id,x=20,y=0,z=0))
            
            
            
            net.populations.append(exc_pop)
            net.populations.append(inh_pop)
            
            exc_syn = Synapse(id='rsExc', lems_source_file='figure1c_Parameters.xml')
            inh_syn = Synapse(id='rsInh', lems_source_file='figure1c_Parameters.xml')
            net.synapses.append(exc_syn)
            net.synapses.append(inh_syn)
            
            #### Weak GBA
            syns = {exc_pop.id:exc_syn.id, inh_pop.id:inh_syn.id}
            
            W = [['wee', 'wie'],
                 ['wei','wii']]
            
            # Add internal connections
            pops = [exc_pop, inh_pop]
            internal_connections(pops,W,syns)
            
            
            # Save to JSON format
            net.id = 'Joglekar_figure1b'
            new_file = net.to_json_file('Joglekar_figure1c.json')
            
            sim = Simulation(id='SimJoglekar_figure1c',
                                                duration='2',
                                                dt='0.02',
                                                network=new_file,
                                                record_variables={'r':{'all':'*'}}
                                                )
                                        
            sim.to_json_file('SimJoglekar_figure1c.nmllite.json')
            
            check_to_generate_or_run(sys.argv,sim)
            
            # Open data file 
            data = np.loadtxt('Excitatory_'+str(round(wEE,5))+'_'+str(round(wEI,5))+'_0.r.dat')
            maxamplocal[j,i] = np.max(data[:,1])
            os.remove('Excitatory_'+str(round(wEE,5))+'_'+str(round(wEI,5))+'_0.r.dat')
            os.remove('Inhibitory_'+str(round(wEE,5))+'_'+str(round(wEI,5))+'_0.r.dat')
            
        else:
            # set value -5 to unstable region (it helps to plot the colormap) 
            maxamplocal[j,i]=-5
            
            
np.save('maxamplocal.npy',maxamplocal)           
