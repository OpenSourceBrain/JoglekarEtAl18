import sys

sim = 'neuron'
print(sys.argv)
if len(sys.argv)>1: 
    sim = sys.argv[1]
    sys.argv = [sys.argv[0]]

sys.argv.append('2')
sys.argv.append('5B')
sys.argv.append('-%s'%sim)
sys.argv.append('yes')
sys.argv.append('yes')
sys.argv.append('no')


import figuresPyNN