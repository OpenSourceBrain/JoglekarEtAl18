import sys

if '-nogui' in sys.argv:
    sys.argv.remove('-nogui')

sys.argv.append('2')
sys.argv.append('5B')
sys.argv.append('yes')
sys.argv.append('yes')
sys.argv.append('no')


import figures