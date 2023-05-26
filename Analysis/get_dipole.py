# Example usage:
# python get_dipole.py -in trajectories/plus100ns_start_water_300K_1to025_id3_r0_q1_rep0_4_output.h5 -out test -l L13

import mdtraj
import numpy as np
import tqdm
import argparse
import os
from openmm.app import AmberPrmtopFile
from scipy.spatial.distance import cdist

parser = argparse.ArgumentParser(description='Run ML MM')
parser.add_argument('-in','--infile',type=str,help='Infile')
parser.add_argument('-out','--outfile',type=str,help='Outfile')
parser.add_argument('-l','--linker',type=str,help='Linkername')
args = parser.parse_args()

assert os.path.isfile(args.infile)
assert not os.path.isfile(args.outfile)

pairdic = {'L2':(10,259),'L6':(24,275),'L9':(11,278),'L13':(11,294)}
kappadic = {'L2':(10,13),'L6':(24,21),'L9':(11,14),'L13':(11,14)}


def get_distances(xyz,linkername):
    disid = pairdic[linkername]
    return np.sqrt(np.sum((xyz[:,disid[0]] - xyz[:,disid[1]]) ** 2,axis=1))


dye_list = {'L2':'A_Dye1_plus_helix_cu','L6':'K_Dye1_plus_helix_cu','L9':'A_Dye2_plus_helix_cu','L13':'K_Dye2_plus_helix_cu'}

top_file = '../Simulation/topologies_and_starting_coordinates/'
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

distances = np.array(distances)
kappa = np.array(kappa)
np.save(args.outfile + '_dis.npy', distances)
np.save(args.outfile + '_kappa.npy',kappa)
np.save(args.outfile + '_dipole.npy', dipoles)
