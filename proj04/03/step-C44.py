
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

name = 'C44'

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

OPTIONS = np.linspace(0.00, 0.05, 6)
volumes = []
energies = []

cr = 0.18
mn = 0.10
fe = 1.0-cr-mn

for opt in OPTIONS:

    l = 3.597/np.sqrt(2)
    atoms = bulk('Fe', 'bcc', a=l)

    scale = [[1,0,0],[0,1,0],[0,0,np.sqrt(2)]]
    atoms.set_cell(np.dot( atoms.get_cell(), scale), scale_atoms=True)

    atoms.set_tags([1])

    alloys = []
    alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))
    alloys.append(Alloy(1, 'Cr', cr / 2, 1.0))
    alloys.append(Alloy(1, 'Cr', cr / 2, -1.0))
    alloys.append(Alloy(1, 'Mn', mn / 2, 1.0))
    alloys.append(Alloy(1, 'Mn', mn / 2, -1.0))

#    dist = [[1+opt, 0, 0], [0, 1+opt, 0], [0, 0, 1/(1+opt)**2]]
    dist = [[1+opt, 0 , 0], [0, 1-opt, 0], [0, 0, 1 / (1 - opt ** 2)]]

    atoms.set_cell(np.dot( atoms.get_cell(), dist), scale_atoms=True)

    calc = EMTO()
    calc.set(dir='{0}/calc/{1}/opt-{2:0.3f}'.format(temp_dir, name, opt),
             lat=10,
             kpts=[13, 13, 13],
             dmax=2.40
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


OPTIONS = (-1.0*OPTIONS[::-1]).tolist() + OPTIONS.tolist()
energies = (energies[::-1]) + energies


coefs = poly.polyfit(OPTIONS, energies, 3)

C44 = coefs[2]/volumes[0]/kJ*1.0e24

print C44

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

save(result_sum, '{0}, {1}, {2}, {3}'.format(name, C44, volumes, energies))





