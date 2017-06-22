
import csv
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk
from mpi4py import MPI
from ase.units import *
import numpy.polynomial.polynomial as poly

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print rank, size

name = '16'

curr_dir = os.getcwd()

temp_dir = curr_dir.replace('work','temp')

os.system('mkdir -p {0}'.format(temp_dir))

result = '{0}/result-{1}.txt'.format(temp_dir,name)
os.system('rm {0}'.format(result))

result_sum = '{0}/summary-{1}.csv'.format(temp_dir,name)
os.system('rm {0}'.format(result_sum))

save(result, '{0}'.format(name))
save(result_sum, '{0}'.format(name))

# OPTIONS = np.linspace(0.98, 1.02, 9)
OPTIONS = range(3, 23, 2)
volumes = []
energies = []

print OPTIONS

cr = 0.15
fe = 1.0-cr

alloys = []
alloys.append(Alloy(1, 'Fe', fe / 2.0 , 1.0))
alloys.append(Alloy(1, 'Fe', fe / 2.0 , -1.0))
alloys.append(Alloy(1, 'Cr', cr / 2.0 , 1.0))
alloys.append(Alloy(1, 'Cr', cr / 2.0 , -1.0))

for opt in OPTIONS:
    l = 3.60
    a0 = l / np.sqrt(2)
    c0 = np.sqrt(8 / 3.0) * a0

    atoms = Atoms('Fe3',
                  scaled_positions=[(0, 0, 0),
                                    (1. / 3., 1. / 3., 1. / 3.),
                                    (2. / 3., 2. / 3., 2. / 3.)],
                  cell=[[1. / 2., sqrt(3) / 2., 0], [-1. / 2., sqrt(3) / 2., 0], [0, 0, 1.5]],
                  pbc=(1, 1, 1))

    scale = [[a0, 0, 0], [0, a0, 0], [0, 0, c0]]
    atoms.set_cell(np.dot(atoms.get_cell(), scale), scale_atoms=True)

    atoms.set_tags([1, 1, 1])

    calc = EMTO()
    calc.set(dir='{0}/{1}/opt-{2}'.format(temp_dir, name, opt),
             lat=9,
             kpts=[opt,opt,opt],
             amix=0.02,
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)

    nm_e = 0
    nm_e = atoms.get_potential_energy() / atoms.get_number_of_atoms()
    nm_v = atoms.get_volume() / atoms.get_number_of_atoms()

    if nm_e < -0.001:
        volumes.append(nm_v)
        energies.append(nm_e)

    save(result, '{3} result : {0} {1} {2}'.format(opt, nm_v, nm_e, name))


save(result, OPTIONS)
save(result, energies)

save(result, '------------------------')

save(result_sum, '{0}, {1}, {2}'.format(name, volumes, energies))

plt.clf()
plt.plot(OPTIONS, energies)
plt.xlabel('OPTIONS')
plt.ylabel('Energy (eV/atom)')

plt.savefig('step-{0}.png'.format(id))

os.system('mv step-{0}.png {1}'.format(id, curr_dir))









