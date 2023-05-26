'''
File to setup Replica exchange run
'''

from sys import stdout
import openmm.app as app
import openmm.unit as unit
import openmm as mm
from openmm.app import NoCutoff, HBonds
from Forces.Nonbonded import *
from mdtraj.reporters import HDF5Reporter
from copy import deepcopy
import numpy as np
from openmm.vec3 import Vec3
from openmm.unit import nanometer, Quantity


class ReplicaExchange:

    def __init__(self, prmtop, inpcrd, number_of_replicas=1):

        self._prmtop = prmtop
        self._inpcrd = inpcrd
        self._numrep = number_of_replicas

        # Create simulation
        self._system = prmtop.createSystem(nonbondedMethod=NoCutoff,nonbondedCutoff=None, constraints=HBonds)

        # List of states
        self._states = [None for i in range(self._numrep)]

        # List of Parameters
        self._parameters = [None for i in range(self._numrep)]

        # List of Energies
        self._energies = [None for i in range(self._numrep)]

        # List of velocities
        self._velocities = [None for i in range(self._numrep)]

        self._dof = self.calculate_dof()
        self._replica_exchange_changes = []
        self._runname = 'dummy_name'

    def calculate_dof(self):
        dof = 0
        for i in range(self._system.getNumParticles()):
            if self._system.getParticleMass(i) > 0 * unit.dalton:
                dof += 3
        for i in range(self._system.getNumConstraints()):
            p1, p2, distance = self._system.getConstraintParameters(i)
            if self._system.getParticleMass(p1) > 0 * unit.dalton or self._system.getParticleMass(p2) > 0 * unit.dalton:
                dof -= 1
        if any(type(self._system.getForce(i)) == mm.CMMotionRemover for i in range(self._system.getNumForces())):
            dof -= 3
        return dof

    def replace_nonbonded_force_with_custom(self):

        for custom_force in [Custom_electrostatic, Custom_lennard_jones, Custom_exception_force_with_scale]:
            cf = custom_force()
            cf.get_particles_from_existing_nonbonded_force(self._system.getForces()[3])
            self._system.addForce(cf._force)

        for i, m in enumerate(self._system.getForces()):
            if isinstance(m, mm.NonbondedForce):
                self._system.removeForce(i)
                break
    
    def create_opencl_simulation(self,positions=None):
        integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        platform = mm.Platform.getPlatformByName('OpenCL')
        simulation = app.Simulation(self._prmtop.topology, self._system, integrator, platform)
        if positions is None:
            simulation.context.setPositions(self._inpcrd.positions)
        else:
            simulation.context.setPositions(positions)
        self._simulation = simulation

    def create_cpu_simulation(self,positions=None):
        integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        platform = mm.Platform.getPlatformByName('CPU')
        simulation = app.Simulation(self._prmtop.topology, self._system, integrator, platform)
        if positions is None:
            simulation.context.setPositions(self._inpcrd.positions)
        else:
            simulation.context.setPositions(positions)
        self._simulation = simulation

    def create_simulation(self,temperature=300,positions=None):

        integrator = mm.LangevinMiddleIntegrator(temperature * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        platform = mm.Platform.getPlatformByName('CUDA')
        platformProperties = {'Precision': 'mixed', 'CudaPrecision': 'mixed'}
        simulation = app.Simulation(self._prmtop.topology, self._system, integrator, platform,platformProperties)
        if positions is None:
            simulation.context.setPositions(self._inpcrd.positions)
        else:
            simulation.context.setPositions(positions)
        self._simulation = simulation

    def add_reporters(self,n_interval = 1000,runname='test',iteration=0):
        self._runname = runname
        self._simulation.reporters.append(app.StateDataReporter('trajectories/' + runname + '_'+ str(iteration) +  '_info.txt', n_interval, step=True,
                                                         potentialEnergy=True, kineticEnergy=True, temperature=True,
                                                         speed=True))
        self._simulation.reporters.append(HDF5Reporter('trajectories/' +runname + '_' + str(iteration)+ '_output.h5', n_interval))

    def minimize(self):

        self._simulation.minimizeEnergy()

    def get_state(self):

        return self._simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True,
                                                  getParameters=True, getParameterDerivatives=True,
                                                  getIntegratorParameters=True)

    def save_states(self,iteration=0):

        for i in range(self._numrep):
            self._simulation.context.setState(self._states[i])
            self._simulation.saveState('states/' + self._runname + '_state_' + str(i) + '_' + str(iteration) + '.xml')
        self._iteration = iteration

    def load_states(self,iteration=0):

        if iteration != 0:
            for i in range(self._numrep):
                self._simulation.loadState('states/' + self._runname + '_state_' + str(i) + '_' + str(iteration) + '.xml')
                self._states[i] = self.get_state()
        else:
            print('state were randomly initialized')

    def initialize_states(self):
        
        state = self.get_state()
        self._states = [deepcopy(state) for i in range(self._numrep)]
        self._velocities = self.get_velocities_of_states()

    def set_positions_based_on_parameters(self,dir,ran,id,q):

        for i in range(self._numrep):

            if i > 10:
                load_id = 10
            else:
                load_id = i
            pos = np.load(dir + 'rep%i_ran%i_id%i_qa%i.npy' % (load_id,ran,id,q))
            pos_q = to_Vec3_quant(pos[0])
            print(pos_q)
            self._simulation.context.setState(self._states[i])
            self._simulation.context.setPositions(pos_q)
            self._states[i] = self.get_state()

    def set_parameters(self,parameters):

        self._parameters = parameters

    def set_original_states(self):

        self._states = []

    def get_total_energy_of_state(self,state):

        return state.getKineticEnergy() + state.getPotentialEnergy()

    def make_simulation_step(self,nsteps=100):

        for i in range(self._numrep):

            self._simulation.context.setState(self._states[i])
            self.set_parameter_from_idx(i)
            self._simulation.step(nsteps)
            self._states[i] = self.get_state()

    def set_parameter_name(self,name):
        self._param_name = name

    def set_parameter_from_idx(self,idx):

        self._simulation.context.setParameter(self._param_name, self._parameters[idx])

    def get_total_energies_of_states(self):

        return [self.get_total_energy_of_state(state) for state in self._states]

    def get_pot_energy_of_state(self,state):

        return state.getPotentialEnergy()

    def get_pot_energies_of_states(self):

        return [self.get_pot_energy_of_state(state) for state in self._states]

    def get_temperature_of_state(self,state):
        return 2*state.getKineticEnergy()/(self._dof*unit.MOLAR_GAS_CONSTANT_R)

    def get_temperatures_of_states(self):

        return [self.get_temperature_of_state(state) for state in self._states]

    def calculate_exchange_probability(self,ene1,ene2,T1,T2):

        k = unit.BOLTZMANN_CONSTANT_kB
        print(T1,T2)
        print((ene1-ene2),(1/(k*T1) - 1/(k*T2)))

        expo = (ene1-ene2)*(1/(k*T1) - 1/(k*T2))
        expo_per_particle = expo/unit.AVOGADRO_CONSTANT_NA
        print(expo_per_particle)
        exp_frac = np.exp(expo_per_particle)
        p = min(1.0,exp_frac)

        return p

    def get_velocity_of_state(self,state):

        return state.getVelocities()

    def get_velocities_of_states(self):

        return [self.get_velocity_of_state(state) for state in self._states]

    def try_exchange(self,indices):
        '''
        Implementation of the Metropolise criterion
        Formulas taken from:
        https://aip.scitation.org/doi/pdf/10.1063/1.1308516
        :param indices:
        :return:
        '''

        temperatures = self.get_temperatures_of_states()

        for i,j in indices:
            self._simulation.context.setState(self._states[i])
            # analyse position x1
            self.set_parameter_from_idx(i)
            h1x1 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()
            self.set_parameter_from_idx(j)
            h2x1 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()

            # analyse position x2
            self._simulation.context.setState(self._states[j])
            self.set_parameter_from_idx(i)
            h1x2 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()
            self.set_parameter_from_idx(j)
            h2x2 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()

            # calculate beta
            beta1 = 1 / (unit.BOLTZMANN_CONSTANT_kB * temperatures[i])
            beta2 = 1 / (unit.BOLTZMANN_CONSTANT_kB * temperatures[j])

            # calculate probability
            deltaH1 = h1x1 - h1x2
            deltaH2 = h2x1 - h2x2

            p = np.exp((beta1 * deltaH1 - beta2 * deltaH2) / unit.AVOGADRO_CONSTANT_NA)
            prob = min(1,p)
            # Switch if probability larger random number
            if prob > np.random.rand(1):
                switch_state = self._states[i]
                self._states[i] = self._states[j]
                self._states[j] = switch_state
                self._replica_exchange_changes.append((i,j))

    def try_exchange_const_t(self,indices,temp=300):
        '''
        Implementation of the Metropolise criterion
        Formulas taken from:
        https://aip.scitation.org/doi/pdf/10.1063/1.1308516
        :param indices:
        :return:
        '''

        temp = temp * unit.kelvin

        for i,j in indices:
            self._simulation.context.setState(self._states[i])
            # analyse position x1
            self.set_parameter_from_idx(i)
            h1x1 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()
            self.set_parameter_from_idx(j)
            h2x1 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()

            # analyse position x2
            self._simulation.context.setState(self._states[j])
            self.set_parameter_from_idx(i)
            h1x2 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()
            self.set_parameter_from_idx(j)
            h2x2 = self._simulation.context.getState(getEnergy=True).getPotentialEnergy()

            # calculate beta
            beta1 = 1 / (unit.BOLTZMANN_CONSTANT_kB * temp)
            beta2 = 1 / (unit.BOLTZMANN_CONSTANT_kB * temp)

            # calculate probability
            deltaH1 = h1x1 - h1x2
            deltaH2 = h2x1 - h2x2

            p = np.exp((beta1 * deltaH1 - beta2 * deltaH2) / unit.AVOGADRO_CONSTANT_NA)
            prob = min(1,p)
            # Switch if probability larger random number
            if prob > np.random.rand(1):
                switch_state = self._states[i]
                self._states[i] = self._states[j]
                self._states[j] = switch_state
                self._replica_exchange_changes.append((i,j))

    def try_temperature_exchange(self, indices):

        energies = self.get_pot_energies_of_states()
        temperatures = self.get_temperatures_of_states()
        new_indices = [i for i in range(self._numrep)]
        self._velocities = self.get_velocities_of_states()
        new_velocities = deepcopy(self._velocities)
        for i,j in indices:
            # Calculate probability
            p = self.calculate_exchange_probability(energies[i],energies[j],temperatures[i],temperatures[j])
            # check probability against random number
            if p > np.random.rand(1):
                new_indices[i] = j
                new_indices[j] = i
                factor = np.sqrt(temperatures[i] / temperatures[j])
                new_velocities[i] = self._velocities[i] * factor
                new_velocities[j] = self._velocities[j] / factor

        self._velocities = new_velocities
        self._states = [self._states[idx] for idx in new_indices]

    def save_exchanges(self):
        arr = np.array(self._replica_exchange_changes)
        np.save('states/' + self._runname + '_' + str(self._iteration) ,arr)


def to_Vec3(xyz):
    return Vec3(x=xyz[0], y=xyz[1], z=xyz[2])


def to_Vec3_quant(pos):
    vectors = []
    for p in pos:
        vectors += [to_Vec3(p)]
    return Quantity(vectors, nanometer)


def get_values(value):
    return [convert_unit_to_float(value[0]), convert_unit_to_float(value[1]), convert_unit_to_float(value[2])]


def convert_unit_to_float(unit):
    try:
        fl = unit._value
    except:
        fl = unit

    return fl

