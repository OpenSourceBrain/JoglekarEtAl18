from neuromllite.NetworkGenerator import check_to_generate_or_run
from neuromllite import Simulation

from neuromllite import Network, Population, Projection, Cell, Synapse, InputSource, Input
from neuromllite import RandomConnectivity,RectangularRegion, RelativeLayout

import sys
import numpy as np

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
                


                
# Build the network
net = Network(id='Joglekar_figure1b')
net.notes = 'A simple rate model with E and I populations'


r1 = RectangularRegion(id='Joglekar', x=0,y=0,z=0,width=1000,height=100,depth=1000)
net.regions.append(r1)



wEEweak=4.45                                # synpatic weight E to E weak LBA    
wEIweak=4.7                                 # synpatic weight I to E weak LBA
wEEstrong=6                                 # synaptic weight E to E strong LBA
wEIstrong=6.7                               # synaptic weight I to E strong LBA
wIE=4 + 2.0/7.0                             # synpatic weight E to I
wII=wIE*1.1                                 # synpatic weight I to I 
tE=0.02 

Astrong=(1/tE)*np.array([[wEEstrong-1, -wEIstrong],[wIE, -wII-1]])
Aweak=(1/tE)*np.array([[wEEweak-1, -wEIweak],[wIE, -wII-1]])

net.parameters = { 'weeStrong':      Astrong[0,0],
                   'weiStrong':      Astrong[1,0],
                   'wieStrong':      Astrong[0,1],
                   'wiiStrong':      Astrong[1,1],
                   'weeWeak':      Aweak[0,0],
                   'weiWeak':      Aweak[1,0],
                   'wieWeak':      Aweak[0,1],
                   'wiiWeak':      Aweak[1,1]}  


exc_cell = Cell(id='Exc', lems_source_file='figure1b_Parameters.xml')
inh_cell = Cell(id='Inh', lems_source_file='figure1b_Parameters.xml')
net.cells.append(exc_cell)
net.cells.append(inh_cell)

exc_popWeak = Population(id='ExcitatoryWeak', 
                     size=1, 
                     component=exc_cell.id, 
                     properties={'color': '0.8 0 0','radius':10},
                     relative_layout = RelativeLayout(region=r1.id,x=-20,y=0,z=0))

inh_popWeak = Population(id='InhibitoryWeak', 
                     size=1, 
                     component=inh_cell.id, 
                     properties={'color': '0 0 0.8','radius':10},
                     relative_layout = RelativeLayout(region=r1.id,x=20,y=0,z=0))

exc_popStrong = Population(id='ExcitatoryStrong', 
                     size=1, 
                     component=exc_cell.id, 
                     properties={'color': '0.8 0 0','radius':10},
                     relative_layout = RelativeLayout(region=r1.id,x=-20,y=0,z=0))

inh_popStrong = Population(id='InhibitoryStrong', 
                     size=1, 
                     component=inh_cell.id, 
                     properties={'color': '0 0 0.8','radius':10},
                     relative_layout = RelativeLayout(region=r1.id,x=20,y=0,z=0))

net.populations.append(exc_popWeak)
net.populations.append(inh_popWeak)
net.populations.append(exc_popStrong)
net.populations.append(inh_popStrong)

exc_syn = Synapse(id='rsExc', lems_source_file='figure1b_Parameters.xml')
inh_syn = Synapse(id='rsInh', lems_source_file='figure1b_Parameters.xml')
net.synapses.append(exc_syn)
net.synapses.append(inh_syn)

#### Weak GBA
syns = {exc_popWeak.id:exc_syn.id, inh_popWeak.id:inh_syn.id}

W = [['weeWeak', 'wieWeak'],
     ['weiWeak','wiiWeak']]

# Add internal connections
pops = [exc_popWeak, inh_popWeak]
internal_connections(pops,W,syns)

#### Strong GBA
syns = {exc_popStrong.id:exc_syn.id, inh_popStrong.id:inh_syn.id}

W = [['weeStrong', 'wieStrong'],
     ['weiStrong','wiiStrong']]

# Add internal connections
pops = [exc_popStrong, inh_popStrong]
internal_connections(pops,W,syns)


# Save to JSON format
net.id = 'Joglekar_figure1b'
new_file = net.to_json_file('Joglekar_figure1b.json')

sim = Simulation(id='SimJoglekar_figure1b',
                                    duration='600',
                                    dt='0.01',
                                    network=new_file,
                                    recordVariables={'r':{'all':'*'}}
                                    )
                            
sim.to_json_file('SimJoglekar_figure1b.nmllite.json')

check_to_generate_or_run(sys.argv,sim)
