"""
November 15, 2018

This script generates an MD simulation. It parses command line arguments for the following parameters and 
produces erg and rvc output files. An input file name is required.

'--iF' [rvc format input file] '--kB' [bonded constant] '--kN' [non-bonded constant] '--nbCutoff' [distance cutoff 
for non-bonded pairs] '--m' [atom mass] '--dt' [time step] '--n' [number of iterations] '--out' [output file stem]

Usage: python3 MyMD.py --iF input.rvc <other params>

"""
import argparse
import numpy as np

class MDSim(object):
    def __init__(self, **kwargs):
        """
        Mothership class for the MD simulation. 
        """
        self.num_atoms = 0
        self.positions = None
        self.velocities = None
        self.accelerations = None
        self.bonded_atoms = {}
        self.nonbonded_atoms = {}
        self.__dict__.update(kwargs) # Store parsed parameters as class attributes
        self.ref_dist_matrix = None
        self.dist_matrix = None
        self.initial_energy = 0
        self.PEb = 0
        self.PEn = 0
        self.KE = 0
        self.force_array = None
        self.vel_t_dt2 = None
        self.inputLines = []
        self.total_energy = 0

    def parse_rvc_file(self):
        """
        Parses an RVC (positions/velocities/connectivities) file.

        Parameters: 
        -----------
        None

        Returns: 
        --------
        None
        """
        # Read in the input rvc file
        input_file = self.iF
        f = open(input_file, 'r')

        self.inputLines = [line.strip() for line in f if line.strip()]
        positions_list = []
        velocities_list = []
        for line in self.inputLines[1:]:
            atom_ID, x_pos, y_pos, z_pos, x_vel, y_vel, z_vel = line.split()[:7]
            positions_list.append([x_pos, y_pos, z_pos])
            velocities_list.append([x_vel, y_vel, z_vel])
            bonded_atoms = list(int(i) for i in line.split()[7:])
            # Store bonded atoms as a dictionary with key, value as atom, list(bonded atoms)
            self.bonded_atoms[int(atom_ID)] = bonded_atoms

        f.close()

        # Assign attributes; convert all lists to numpy arrays
        self.num_atoms = len(self.bonded_atoms)
        self.accelerations = np.zeros( (self.num_atoms, 3) )
        self.positions = np.array(positions_list, dtype='float_')
        self.velocities = np.array(velocities_list, dtype='float_')

        # Initialize the system energy with KE (PE=0 at t=0)
        v_squared = np.sum( self.velocities**2, axis=1 )
        self.initial_energy = 0.5 * self.m * np.sum((v_squared))

        # Compute Euclidean distances between all atoms to identify fixed set of non-bonded atoms and 
        # all reference distances. Fill ref_dist_matrix with pairwise distances (b(0)). dist_matrix
        # will be populated with distances (b) after updating the positions.
        self.ref_dist_matrix = np.zeros( (self.num_atoms, self.num_atoms) )
        self.dist_matrix = np.zeros( (self.num_atoms, self.num_atoms) )

        for i in range( self.num_atoms ):
            # Fill the diagonal with arbitrary large number to avoid non-bonded self-pair assignment
            self.ref_dist_matrix[i][i] = 10.0**12
            for j in range(i+1, self.num_atoms):
                first_array = self.positions[i,:]
                second_array = self.positions[j,:]
                distance = np.linalg.norm(first_array-second_array)
                self.ref_dist_matrix[i][j] = distance
                self.ref_dist_matrix[j][i] = distance

        # Store non-bonded atoms as a dictionary with key, value as atom, list(non-bonded atoms)
        # Non-bonded atoms are atoms less than nbCutoff that are not also bonded to a given atom
        for i in range(self.num_atoms):
            nonbonded_atom_index = np.where(self.ref_dist_matrix[i] < self.nbCutoff)[0]
            self.nonbonded_atoms[i+1] = [x+1 for x in nonbonded_atom_index if (not x+1 in self.bonded_atoms[i+1]) and (x != i)]

    def do_verlet_iteration(self):
        """
        Performs one iteration of the Velocity Verlet integrator.
        Inputs:
            None
        Returns:
            None
        """
        self.vel_t_dt2 = self.velocities + 0.5 * self.accelerations * self.dt
        self.update_positions()
        self.compute_force()
        self.accelerations = (1./self.m) * self.force_array
        self.velocities = self.vel_t_dt2 + 0.5 * self.accelerations * self.dt # not right

        # Compute kinetic energy
        v_squared = np.sum( self.velocities**2, axis=1 )
        self.KE = 0.5 * self.m * np.sum((v_squared))

        # Compute total energy
        self.total_energy = self.PEb + self.PEn + self.KE

    def update_positions(self):
        """
        Updates atom positions according to the Velocity Verlet integrator. 
        Also updates distance matrix (pairwise Euclidean distances) using new positions.

        Inputs:
            None
        Returns:
            None
        """
        self.positions += (self.vel_t_dt2 * self.dt)

        for atom_A in range(1, self.num_atoms+1):
            for atom_B in self.bonded_atoms[atom_A]:
                distance = np.linalg.norm(self.positions[atom_A-1,:] - self.positions[atom_B-1,:])
                self.dist_matrix[atom_A-1][atom_B-1] = distance
                self.dist_matrix[atom_B-1][atom_A-1] = distance

            for atom_C in self.nonbonded_atoms[atom_A]:
                distance = np.linalg.norm(self.positions[atom_A-1,:] - self.positions[atom_C-1,:])
                self.dist_matrix[atom_A-1][atom_C-1] = distance
                self.dist_matrix[atom_C-1][atom_A-1] = distance

    def compute_force(self):
        '''
        Computes the force on a given atom from every other interacting atom (bonded and non-bonded).
        Also computes potential energy for the system.
        Inputs:
            None
        Returns:
            None
        '''
        self.PEb = 0
        self.PEn = 0
        # Initialize an array to hold forces in x, y, z on each atom from every other interacting atom
        self.force_array = np.zeros( (self.num_atoms, 3) )

        for atom_A in range(1, self.num_atoms+1):
            for atom_B in self.bonded_atoms[atom_A]:

                distance = self.dist_matrix[atom_A-1][atom_B-1] - self.ref_dist_matrix[atom_A-1][atom_B-1]
                magnitude = self.kB * distance
                displacement = self.positions[atom_B-1,:] - self.positions[atom_A-1,:]
                force = displacement * magnitude * 1.0/(self.dist_matrix[atom_A-1][atom_B-1])
                self.force_array[atom_A-1,:] += force
                
                self.PEb += (0.5 * self.kB * distance**2)

            for atom_C in self.nonbonded_atoms[atom_A]:

                distance = self.dist_matrix[atom_A-1][atom_C-1] - self.ref_dist_matrix[atom_A-1][atom_C-1]
                magnitude = self.kN * distance
                displacement = self.positions[atom_C-1,:] - self.positions[atom_A-1,:]
                force = displacement * magnitude * 1.0/(self.dist_matrix[atom_A-1][atom_C-1])
                self.force_array[atom_A-1,:] += force
                
                self.PEn += (0.5 * self.kN * distance**2)
                
        # Correct potential energies for reciprocal interactions
        self.PEb /= 2.0
        self.PEn /= 2.0

    def write_rvc_output(self, iter_num, file_handle):
        """
        Writes current atom positions and velocities, as well as the iteration
        number, to an output file given by file_handle.

        Parameters: 
        -----------
        iter_num : int, current iteration number.
        file_handle : handle to *_out.rvc file opened for writing
                      (*NOT* file name)

        Returns:
        --------
        None        
        """
        # Copy the first frame from the input rvc file
        if iter_num == 0:
            file_handle.write(self.inputLines[0])
            for line in self.inputLines[1:]:
                file_handle.write('\n')
                file_handle.write('\t'.join(list(map(str,line.split()))))
        # Otherwise, add new frame
        else:
            file_handle.write('\n# At time step ' + str(iter_num) +', energy = ' + '{:.3f}'.format(self.total_energy) + 'kJ')
            for i in range(self.num_atoms):
                file_handle.write('\n'+str(i+1)+'\t')
                out_list = ['{:.4f}'.format(k) for k in self.positions[i,:].tolist()] + ['{:.4f}'.format(j) for j in self.velocities[i,:].tolist()] + [str(g) for g in self.bonded_atoms[i+1]]
                file_handle.write('\t'.join( out_list))
        
    def write_erg_output(self, iter_num, file_handle):
        """
        Writes energy statistics (kinetic energy, potential energy of 
        bonded interactions, potential energy of nonbonded interactions,
        and the sum of the foregoing energies - E_tot) as well as the iteration
        number to an output file given by file_handle.

        Parameters: 
        -----------
        iter_num : int, current iteration number.
        file_handle : handle to *_out.erg file opened for writing 
                     (*NOT* file name)

        Returns:
        --------
        None
        """
        out_list = ['{:.1f}'.format(k) for k in [self.KE, self.PEb, self.PEn, self.total_energy]]
        file_handle.write('\t'.join( ['\n' + str(iter_num)] + out_list ))

    def run_md_sim(self):
        """
        Runs the MD simulation.
        Inputs:
            None
        Returns:
            None
        """
        rvc_file = open(self.out + '_out.rvc', 'w')
        erg_file = open(self.out + '_out.erg', 'w')

        # Add header line to erg file
        erg_file.write('\t'.join( ['# step', 'E_k', 'E_b', 'E_nB', 'E_tot'] ))
        # Copy the input rvc file as the first frame of the rvc output file
        self.write_rvc_output(0, rvc_file)

        # Run the simulation for n iterations
        for i in range(1, self.n+1):
            self.do_verlet_iteration()

            # Check for overflow errors (fluctuating energy)
            # Print an error statement and terminate program if found
            if (self.initial_energy * 10) < self.total_energy or self.total_energy < (self.initial_energy/10):
                print('Overflow error at step ' + str(i) +'. Simulation terminated.')
                self.write_rvc_output(i, rvc_file)
                self.write_erg_output(i, erg_file)
                rvc_file.close()
                erg_file.close()
                return

            # On every 10th iteration, print to the output files
            if i % 10 == 0:
                self.write_rvc_output(i, rvc_file)
                self.write_erg_output(i, erg_file)

        rvc_file.close()
        erg_file.close()

class ParsedArgs(object):
    '''
    Class to hold the argparse parser, apply default values for arguments, and add arguments to
    the parser from the command line.
    '''
    def __init__(self):
        # Initialize the parser object and default values
        self.parser = argparse.ArgumentParser()
        parameter_list = ['--iF', '--kB', '--kN', '--nbCutoff', '--m', '--dt', '--n']
        parameter_default = [None, 40000.0, 400.0, 0.50, 12.0, 0.001, 1000]
        parameter_type = [str, float, float, float, float, float, int]

		# Add arguments to parser (except --out, which depends on first parsing --iF)
        # If --out was provided it will not be in 'namespace' until added as an argument
        for i, j, k in zip(parameter_list, parameter_default, parameter_type):
            self.fill_parser(i, j, k)
        namespace, extra = self.parser.parse_known_args()
		
		# Make sure a filename was specified; if not, print an error message and terminate
        if namespace.iF == None:
            print("Error: must specify an input file")
            return

		# Add the --out argument to the parser
        self.fill_parser('--out', namespace.iF.split('.')[0], str)
        self.parsed_parameters = vars(self.parser.parse_args())

    def fill_parser(self, parameter, p_default, p_type):
        '''
        Adds arguments to the parser object.
        Inputs:
            parameter = str, parameter flag
            p_default = nonetype or int, default value for parameter
            p_type = type, Python type to assign to value from command line
        Returns:
            None
        '''
        self.parser.add_argument(parameter, action='store', dest=parameter.split('-')[-1], default=p_default, type=p_type)

def main():
    '''
    Runs the script!
    '''
    parameters = ParsedArgs()
    MyFirstMDSim = MDSim(**parameters.parsed_parameters)
    MyFirstMDSim.parse_rvc_file()
    MyFirstMDSim.run_md_sim()


if __name__ == '__main__':
	main()

