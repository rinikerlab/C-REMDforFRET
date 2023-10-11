import mdtraj
import numpy as np
import tqdm
import argparse
import os
from openmm.app import AmberPrmtopFile
from scipy.spatial.distance import cdist

parser = argparse.ArgumentParser(description='Run Analysis')
parser.add_argument('-in','--infile',type=str,help='Infile')
parser.add_argument('-out','--outfile',type=str,help='Outfile')
parser.add_argument('-l','--linker',type=str,help='Linkername')
args = parser.parse_args()

assert os.path.isfile(args.infile)
#assert '.npy' in args.outfile
assert not os.path.isfile(args.outfile)

kappadic = {'0':(340,344),'1':(341,345),'2':(342,346)}
dye_list = {'0':'totdye_hisn','1':'totdye_hisn_3p','2':'totdye_hisn_4p'}

top_file = '/cluster/home/kpaul/whome/projects/lrep/beta/'
top_file += dye_list[args.linker] + '.top'
top = AmberPrmtopFile(top_file)
system = top.createSystem()
nonbonded_force = system.getForces()[3]
charges = np.array([nonbonded_force.getParticleParameters(i)[0]._value for i in range(nonbonded_force.getNumParticles())])

def get_angle(vector_1,vector_2):
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(dot_product)
    return angle

def calculate_kappa_for_frame(xyz,charges,c6_idx,c10_idx):
    copper_dis = cdist(xyz,np.expand_dims(xyz[-1,:],0)).T[0]
    copper_vec = xyz - np.expand_dims(xyz[-1,:],0)
    charges = charges[:-1]
    copper_dis = copper_dis[:-1]
    copper_vec = copper_vec[:-1,:]
    dipole = copper_vec * np.expand_dims(charges / (copper_dis**3),1)
    induced_dipol = np.sum(dipole,0) * -1
    dye_dipole = xyz[c6_idx] - xyz[c10_idx]
    copperC6vec = xyz[c6_idx] - xyz[-1]
    theta = get_angle(induced_dipol,dye_dipole)
    alpha = get_angle(copperC6vec,induced_dipol)
    delta = get_angle(copperC6vec,dye_dipole)
    kappa = (np.cos(theta) - 3*np.cos(delta)*np.cos(alpha))**2
    return kappa , copper_dis[c6_idx], induced_dipol


distances = []
kappa = []
dipoles = []
for i,chunk in tqdm.tqdm(enumerate(mdtraj.iterload(filename=args.infile,chunk=10000))):
    for xyz in chunk.xyz:
        kap, dis, dipole = calculate_kappa_for_frame(xyz,charges,kappadic[args.linker][0],kappadic[args.linker][1])
        distances.append(dis)
        kappa.append(kap)
        dipoles.append(dipole)
#print('trajectory broken at it %i' % i)

distances = np.array(distances)
kappa = np.array(kappa)
np.save(args.outfile + '_dis.npy', distances)
np.save(args.outfile + '_kappa.npy',kappa)
np.save(args.outfile + '_dipole.npy', dipoles)
