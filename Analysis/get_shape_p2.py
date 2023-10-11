import mdtraj
import numpy as np
import tqdm
import argparse
import os
from openmm.app import AmberPrmtopFile
from scipy.spatial.distance import cdist
import pandas as pd

parser = argparse.ArgumentParser(description='Run ML MM')
parser.add_argument('-in','--infile',type=str,help='Infile')
parser.add_argument('-out','--outfile',type=str,help='Outfile')
#parser.add_argument('-l','--linker',type=str,help='Linkername')
args = parser.parse_args()

assert os.path.isfile(args.infile)
assert not os.path.isfile(args.outfile)

kappadic = {'0':(340,344),'1':(341,345),'2':(342,346)}

def calculate_impacts(filenames,chunksize=160):
    folder = os.getenv('TMPDIR') + "/" + args.outfile
    impact_executable = '/cluster/project/igc/kpaul/impact'
    data1 = []
    data2 = []
    data3 = []
    rep_step = 10
    reprange = [i for i in range(16)]
    for filename in filenames:
        for it,traj in tqdm.tqdm(enumerate(mdtraj.iterload(filename,chunk=chunksize))):
            if it % 1 != 0:
                continue
            for rep in reprange:
                t = traj[rep*rep_step + rep_step-1] # get last fram before potential exchange
                t.save_pdb(folder+'rep%s.pdb' % rep)
                
            pdbfiles = [folder+'rep%s.pdb ' % rep for rep in reprange]
            substring = '%s -nThreads 8 -rProbe 1.73 -param param.txt -H ' % impact_executable
            for file in pdbfiles:
                substring += file
            substring += '| grep CCS > %sintermediat.txt' % folder
            os.system(substring)
            data1.append(pd.read_table(folder+'intermediat.txt',delim_whitespace=True,skiprows=1,header=None)[3].values)
            data2.append(pd.read_table(folder+'intermediat.txt',delim_whitespace=True,skiprows=1,header=None).values[:,6])
            data3.append(pd.read_table(folder+'intermediat.txt',delim_whitespace=True,skiprows=1,header=None).values[:,7])

    return np.concatenate(data1), np.concatenate(data2), np.concatenate(data3)

d1,d2,d3 = calculate_impacts([args.infile])

np.save(args.outfile + '_data1.npy', d1)
np.save(args.outfile + '_data2.npy', d2)
np.save(args.outfile + '_data3.npy', d3)


