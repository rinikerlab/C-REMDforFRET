# Example usage:
# python run_at_300K.py -n 0.1 -d 0 -f topologies_and_starting_coordinates/ -e 10000 -r plus100ns_start_water_300K_1to025_id0_r0_q1_rep0 -i 0 -ran 0 -qa 1

import tqdm
import sys
sys.path.append('..')

from ReplicaExchange.ReplicaExchange import ReplicaExchange
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, PDBFile
import argparse

# setup parse
parser = argparse.ArgumentParser(description='Run Replica Exchange')
parser.add_argument('-r','--runname',default='dyename',type=str,help='Name to use for saving files')
parser.add_argument('-n','--ns',default=200,type=float,help='Time of simulation in ns')
parser.add_argument('-d','--dyeid',default=0,type=int,help='ID of the Dye to simulate (starting from 0)')
parser.add_argument('-f','--fileloc',default='.',type=str,help='location of input files')
parser.add_argument('-v','--verbose',default=False,type=bool,help='whether to print progress (not recommended on LSF systems)')
parser.add_argument('-e','--exchangefrequency',default=100,type=int,help='The frequency of exchange trials in steps')
parser.add_argument('-i','--iteration',default=0,type=int,help='Iteration of the process for running in shorter queues')
parser.add_argument('-ran','--random',default=0,type=int,help='random')

args = parser.parse_args()

dye_list = ['A_Dye1_plus_helix_cu','K_Dye1_plus_helix_cu','A_Dye2_plus_helix_cu','K_Dye2_plus_helix_cu']
dye_name = dye_list[args.dyeid]

prmtop = AmberPrmtopFile(args.fileloc + dye_name + '.top')
if dye_name == 'K_Dye2_plus_helix_cu':
    from openmm.app import PDBFile
    inpcrd = PDBFile(args.fileloc + dye_name + '.pdb')
else:
    inpcrd = AmberInpcrdFile(args.fileloc + dye_name + '.crd')


if args.runname == 'dyename':
    runname = dye_name + '_ef_' + str(args.exchangefrequency) + '_cyc_const_T' + '_open_nocutoff'
else:
    runname = args.runname

re = ReplicaExchange(prmtop,inpcrd,16)
re.replace_nonbonded_force_with_custom()
re.create_simulation(300)
re.minimize()

report_frequency = 1000

re.add_reporters(report_frequency,runname,args.iteration)
re.initialize_states()
re.load_states(args.iteration)
epsilon_scale = [1 - i/20 for i in range(16)]
re.set_parameters(epsilon_scale)
re.set_parameter_name('epsilon_factor')

even = [(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15)]
odd = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14)]
index_lists = [even,odd]

steps_in_100 = int(args.ns / 0.002 * 1000 / args.exchangefrequency)

if args.verbose:
    for i in tqdm.tqdm(range(steps_in_100)):
        re.make_simulation_step(args.exchangefrequency)
        re.try_exchange_const_t(index_lists[i % 2],300)
else:
    for i in range(steps_in_100):
        re.make_simulation_step(args.exchangefrequency)
        re.try_exchange_const_t(index_lists[i % 2],300)

re.save_states(args.iteration + 1)
