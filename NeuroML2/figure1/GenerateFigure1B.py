from neuromllite.NetworkGenerator import check_to_generate_or_run
from neuromllite import Simulation

from neuromllite import Network, Population, Projection, Cell, Synapse, InputSource, Input
from neuromllite import RandomConnectivity

import sys
import numpy as np

# This function generates the overview of the network using neuromllite
def internal_connections(pops):
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
                
                
# Build the network
net = Network(id='Joglekar_figure1b')
net.notes = 'A simple rate model with E and I populations'

wEEweak=4.45                                # synpatic weight E to E weak LBA    
wEIweak=4.7                                 # synpatic weight I to E weak LBA
wEEstrong=6                                 # synaptic weight E to E strong LBA
wEIstrong=6.7                               # synaptic weight I to E strong LBA
wIE=4 + 2.0/7.0                             # synpatic weight E to I
wII=wIE*1.1                                 # synpatic weight I to I 
tE=0.02 

Astrong=(1/tE)*np.array([[wEEstrong-1, -wEIstrong],[wIE, -wII-1]])

net.parameters = { 'wee':      Astrong[0,0],
                   'wei':      Astrong[1,0],
                   'wie':      Astrong[0,1],
                   'wii':      Astrong[1,1]}  


exc_cell = Cell(id='Exc', lems_source_file='figure1b_Parameters.xml')
inh_cell = Cell(id='Inh', lems_source_file='figure1b_Parameters.xml')
net.cells.append(exc_cell)
net.cells.append(inh_cell)

exc_pop = Population(id='Excitatory', 
                     size=1, 
                     component=exc_cell.id, 
                     properties={'color': '0.8 0 0','radius':10})

inh_pop = Population(id='Inhibitory', 
                     size=1, 
                     component=inh_cell.id, 
                     properties={'color': '0 0 0.8','radius':10})

net.populations.append(exc_pop)
net.populations.append(inh_pop)

exc_syn = Synapse(id='rsExc', lems_source_file='figure1b_Parameters.xml')
inh_syn = Synapse(id='rsInh', lems_source_file='figure1b_Parameters.xml')
net.synapses.append(exc_syn)
net.synapses.append(inh_syn)

syns = {exc_pop.id:exc_syn.id, inh_pop.id:inh_syn.id}

W = [['wee', 'wie'],
     ['wei','wii']]

# Add internal connections
pops = [exc_pop, inh_pop]
internal_connections(pops)


# Save to JSON format
net.id = 'Joglekar_figure1b'
new_file = net.to_json_file('Joglekar_figure1b.json')

sim = Simulation(id='SimJoglekar_figure1b',
                                    duration='0.6',
                                    dt='0.0001',
                                    network=new_file,
                                    recordRates={'all':'*'}
                                    )
                            
sim.to_json_file('SimJoglekar_figure1b.nmllite.json')

check_to_generate_or_run(sys.argv,sim)
