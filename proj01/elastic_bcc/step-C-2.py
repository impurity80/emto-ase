
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

name = 'C-2'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work','temp')

os.system('mkdir {0}/graph'.format(temp_dir))
os.system('mkdir {0}/result'.format(temp_dir))

result = '{0}/result/result-{1}.txt'.format(temp_dir,name)
os.system('rm {0}'.format(result))

result_sum = '{0}/result/summary-{1}.csv'.format(temp_dir,name)
os.system('rm {0}'.format(result_sum))

save(result, '{0}'.format(name))
save(result_sum, '{0}'.format(name))

OPTIONS = np.linspace(-0.05, 0.05, 11)
volumes = []
energies = []

cr = 0.15
ni = 0.15
fe = 1.0

for opt in OPTIONS:

    l = 2.8784
    atoms = bulk('Fe', 'bcc', a=l)

    atoms.set_tags([1])

    alloys = []
    alloys.append(Alloy(1, 'Fe', fe , 1.0))

    dist = [[1+opt, 0, 0], [0, 1+opt, 0], [0, 0, 1/(1+opt)**2]]
#    dist = [[1+opt, 0 , 0], [0, 1-opt, 0], [0, 0, 1 / (1 - opt ** 2)]]

    atoms.set_cell(np.dot(dist, atoms.get_cell()), scale_atoms=True)

    print atoms.get_cell()

    calc = EMTO()
    calc.set(dir='{0}/calc/{1}/opt-{2:0.3f}'.format(temp_dir, name, opt),
             lat=10,
             kpts=[27,27,27],
             dmax=2.20,
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

print volumes, energies


coefs = poly.polyfit(OPTIONS, energies, 3)

C = coefs[2]/volumes[0]/kJ*1.0e24

print C

x_new = np.linspace(OPTIONS[0]-0.01, OPTIONS[-1]+0.01, num=len(OPTIONS)*10)

ffit = poly.polyval(x_new, coefs)

plt.scatter(OPTIONS, energies)
plt.plot(x_new, ffit)
plt.savefig('{0}.png'.format(name))

os.system('mv {0}.png {1}/graph'.format(name, temp_dir))

save(result, OPTIONS)
save(result, volumes)
save(result, energies)

save(result, '------------------------')

save(result_sum, '{0}, {1}, {2}, {3}'.format(name, C, volumes, energies))







