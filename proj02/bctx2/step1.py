
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
from ase.units import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print rank, size

name = '1'

curr_dir = os.getcwd()

temp_dir = curr_dir.replace('work','temp')

os.system('mkdir -p {0}/graph'.format(temp_dir))
os.system('mkdir -p {0}/result'.format(temp_dir))

result = '{0}/result/result-{1}.txt'.format(temp_dir,name)
os.system('rm {0}'.format(result))

result_sum = '{0}/result/summary-{1}.csv'.format(temp_dir,name)
os.system('rm {0}'.format(result_sum))

save(result, '{0}'.format(name))
save(result_sum, '{0}'.format(name))


OPTIONS = np.linspace(0.98, 1.02, 9)
volumes = []
energies = []

cr = 0.15
fe = 1.0-cr

for opt in OPTIONS:

    l = 3.602 * opt
    a = l / sqrt(2)
    c = l

    atoms = Atoms('Fe2',
                  scaled_positions=[
                      (0.0, 0.0, 0),
                      (0.5, 0.5, 0.5)],
                  cell=[a, a, c],
                  pbc=(1, 1, 1))

    atoms.set_tags([1, 1])

    alloys = []
    alloys.append(Alloy(1, 'Fe', fe, 0.0))
    alloys.append(Alloy(1, 'Cr', cr , 0.0))

    calc = EMTO()
    calc.set(dir='{0}/calc/{1}/opt-{2:0.3f}'.format(temp_dir, name, opt),
             lat=6,
             afm='F',
             kpts=[13, 13, 13],
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)
    nm_e = atoms.get_potential_energy() / atoms.get_number_of_atoms()
    nm_v = atoms.get_volume() / atoms.get_number_of_atoms()

    if nm_e < -0.001:
        volumes.append(nm_v)
        energies.append(nm_e)

    save(result, '{3} result : {0} {1} {2}'.format(opt, nm_v, nm_e, name))

print volumes, energies

eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
eos.plot('{0}/graph/{1}.png'.format(temp_dir, name))

save(result, '{0} {1} {2} {3}'.format(v0, e0, B/kJ*1.0e24, (4.0 * v0) ** (1.0 / 3.0)))

save(result, OPTIONS)
save(result, volumes)
save(result, energies)

save(result, '------------------------')

save(result_sum, '{0}, {1}, {2}, {3}, {4}, {5}'.format(name, e0, v0, B/kJ*1.0e24, volumes, energies))







