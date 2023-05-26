import openmm as mm
from openmm.app import NoCutoff
import openmm.unit as u

class _custom_force:
    def __init__(self):
        self._force = None

    @property
    def openmm_force(self):
        return self._force


class Custom_electrostatic(_custom_force):

    def __init__(self):
        super().__init__()
        # Setup Energy Term
        energy_term = 'epsilon_factor * charge1 * charge2 * fourpieps / r;'

        # Create Force
        force = mm.CustomNonbondedForce(energy_term)

        # Add Global parameters
        force.addGlobalParameter('fourpieps',138.935456)
        force.addGlobalParameter('epsilon_factor',1)

        # Add per particle Parameters
        force.addPerParticleParameter('charge')
        force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)

        self._force = force

    def add_particles(self,charges):

        # Go through charges and set the parameters
        for charge in charges:
            self._force.addParticle([charge])

    def add_exceptions(self,exceptions):

        # Go through charges and set the parameters
        for exception in exceptions:
            self._force.addExclusion(exception[0],exception[1])

    def get_particles_from_existing_nonbonded_force(self,force):

        # Get charges
        for i in range(force.getNumParticles()):
            charge, _, _ = force.getParticleParameters(i)
            self._force.addParticle([charge])

        # Get exclusions
        for i in range(force.getNumExceptions()):
            k, j,_,_,_ = force.getExceptionParameters(i)
            self._force.addExclusion(k,j)


class Custom_exception_force_with_scale(_custom_force):
    '''
    Force to compensate for exceptions in Custom_lennard_jones and Custom_electrostatic
    '''

    def __init__(self):
        super().__init__()

        energy_expression = '4 * epsilon * (sigmar6^2 - sigmar6) + epsilon_factor * chargeprod * fourpieps * (1/r);'
        energy_expression += 'sigmar6 = (sigma/r)^6;'

        force = mm.CustomBondForce(energy_expression)

        force.addGlobalParameter('fourpieps',138.935456)
        force.addGlobalParameter('epsilon_factor',1)

        force.addPerBondParameter('epsilon')
        force.addPerBondParameter('sigma')
        force.addPerBondParameter('chargeprod')

        self._force = force

    def get_particles_from_existing_nonbonded_force(self,force):

        for i in range(force.getNumExceptions()):
            k, j, chargeprod, sigma, epsilon = force.getExceptionParameters(i)
            self._force.addBond(k, j,[epsilon,sigma,chargeprod])


class Custom_lennard_jones(_custom_force):

    def __init__(self):
        super().__init__()
        # Get energy term
        energy_term = '4 * epsilon * ((sigmacom/r)^12 - (sigmacom/r)^6);'
        energy_term += 'sigmacom = 0.5 * (sigma1+sigma2);'
        energy_term += 'epsilon = sqrt(epsilon1*epsilon2);'

        # Create Force
        force = mm.CustomNonbondedForce(energy_term)

        # Add particle parameters
        force.addPerParticleParameter('epsilon')
        force.addPerParticleParameter('sigma')
        force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)

        force.setUseLongRangeCorrection(True)

        self._force = force

    def get_particles_from_existing_nonbonded_force(self,force):

        # Get charges
        for i in range(force.getNumParticles()):
            _, sigma, epsilon = force.getParticleParameters(i)
            self._force.addParticle([epsilon,sigma])

        # Get exclusions
        for i in range(force.getNumExceptions()):
            k, j,_,_,_ = force.getExceptionParameters(i)
            self._force.addExclusion(k,j)