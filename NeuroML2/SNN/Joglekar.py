from neuromllite import Network, Cell, InputSource, Population, Synapse, RectangularRegion, RandomLayout 
from neuromllite import Projection, RandomConnectivity, Input, Simulation
import scipy
import sys

###############################################################################
############################## Parameters #####################################
###############################################################################
R             = 50.0         # membrane resitance [MOhm]
tauE_m        = 20.0         # excitatory neuron membrane time constant [ms]
tauI_m        = 10.0         # inhibitory neuron membrane time constant [ms]
cE_m          = tauE_m/R     # Capacitance excitatory neurons [nF]
cI_m          = tauI_m/R     # Capacitance excitatory neurons [nF]
tauRef        = 2.0          # refractory time [ms]
Vrest         = -70.0        # resting potential [mV]
Vt            = -50.0        # threshold [mV]
Vreset        = -60.0        # reset [mV]  
tau_syn_e     =  0.1         # time constant for synapse [ms] 

# External input for assynchronous behavior
VextE         = 14.2         # External input in excitatory neurons [mV]
VextI         = 14.7         # External input in excitatory neurons [mV]

# Noise
sigmaV=3.0   # [mV]

# Coefficient (delta to alpha synapse)
coeffE=0.67957046
coeffI=2*coeffE

# timestep
dt=0.1
    
# duration
duration=500

# number of areas
nAreas=2

# Axonal velocity
speed=3.5

# Alpha
alpha=4

# mu
muIE     = .19/4 
muEE     = .0375

# # path to .mat files
# path='/home/ronaldo/github/Joglekar2018/Matlab/'
# #hierarchy values file 
# hierVals = scipy.io.loadmat(path+'hierValspython.mat')
# hierValsnew = hierVals['hierVals'][:]
# hier=hierValsnew/max(hierValsnew)#hierarchy normalized. 
# hier=np.squeeze(hier[:nAreas])

# #fln values file 
# flnMatp = scipy.io.loadmat(path+'efelenMatpython.mat')
# conn=flnMatp['flnMatpython'][:][:] #fln values..Cij is strength from j to i 
# conn=conn[:nAreas,:nAreas]

# distMatp = scipy.io.loadmat(path+'subgraphWiring29.mat')
# distMat=distMatp['wiring'][:][:] #distances between areas values..
# delayMat = distMat/speed
# delay=delayMat[:nAreas,:nAreas]

###################### Build the network ######################################
net = Network(id='Joglekar1Network_PyNN')
net.notes = 'Joglekar network: a network with PyNN cells & inputs'

######################## Cell #################################################
cellE = Cell(id='excitatory', pynn_cell='IF_curr_alpha')
cellI = Cell(id='inhibitory', pynn_cell='IF_curr_alpha')

cellE.parameters = {
            "tau_m":tauE_m, 
            "cm":cE_m, 
            "v_rest":Vrest, 
            "v_reset":Vreset, 
            "v_thresh":Vt, 
            "tau_refrac":tauRef,
            "tau_syn_E":tau_syn_e,
            "tau_syn_I":tau_syn_e,
            "i_offset":(VextE/R)}

cellI.parameters = {
            "tau_m":tauI_m, 
            "cm":cI_m, 
            "v_rest":Vrest, 
            "v_reset":Vreset, 
            "v_thresh":Vt, 
            "tau_refrac":tauRef,
            "tau_syn_E":tau_syn_e,
            "tau_syn_I":tau_syn_e,
            "i_offset":(VextI/R)}

### Append cells to network 
net.cells.append(cellE)
net.cells.append(cellI)

######################### Spatial parameters for network ######################
r1 = RectangularRegion(id='region1', x=0,y=0,z=0,width=1000,height=100,depth=1000)
net.regions.append(r1)

############################# Populations #####################################
pE = Population(id='popE', size=1600, component=cellE.id, properties={'color':'1 0 0'},random_layout = RandomLayout(region=r1.id))
pI = Population(id='popI', size=400, component=cellI.id, properties={'color':'1 0 0'},random_layout = RandomLayout(region=r1.id))


### Append populations to network 
net.populations.append(pE)
net.populations.append(pI)

############################### Inputs (Noise) ################################
# Noise on excitatory neurons
stdNoiseE = (sigmaV/R)*(tauE_m**0.5)/(dt**0.5)

input_sourceE = InputSource(id='noisyCurrentSourceE', 
                           lems_source_file='TestNCS.xml', 
                           parameters={'mean':'0.0nA','stdev':str(stdNoiseE)+'nA'})
# Noise on inhibitory neurons    
stdNoiseI = (sigmaV/R)*(tauI_m**0.5)/(dt**0.5)
input_sourceI = InputSource(id='noisyCurrentSourceI', 
                           lems_source_file='TestNCS.xml', 
                           parameters={'mean':'0.0nA','stdev':str(stdNoiseI)+'nA'})


net.input_sources.append(input_sourceE)
net.input_sources.append(input_sourceI)

net.inputs.append(Input(id='stimE',
                        input_source=input_sourceE.id,
                        population=pE.id,
                        percentage=100))
                        
net.inputs.append(Input(id='stimI',
                        input_source=input_sourceI.id,
                        population=pI.id,
                        percentage=100))

####################### Save network in json file #############################
print(net.to_json())
new_file = net.to_json_file('%s.json'%net.id)


################################################################################
###   Build Simulation object & save as JSON

sim = Simulation(id='SimJoglekar1Network',
                 network=new_file,
                 duration='1000',
                 dt='0.01',
                 recordSpikes={pE.id:'*', pI.id:'*'})
                 
sim.to_json_file()


################################################################################
###   Run in some simulators

from neuromllite.NetworkGenerator import check_to_generate_or_run
import sys

check_to_generate_or_run(sys.argv, sim)

 
