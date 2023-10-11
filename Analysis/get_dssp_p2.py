import mdtraj
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Run ML MM')
parser.add_argument('-in','--infile',type=str,help='Infile')
parser.add_argument('-out','--outfile',type=str,help='Outfile')
args = parser.parse_args()

assert os.path.isfile(args.infile)
assert not os.path.isfile(args.outfile)

traj = mdtraj.load(args.infile)
d1 = mdtraj.compute_dssp(traj,False)[:,6:18]
np.save(args.outfile + '_dssp.npy', d1)



