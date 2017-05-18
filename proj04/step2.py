
import os
from emto import *
from emto_eos import *
from ec_fcc import *

from ase.lattice import bulk

id = '2'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)
os.system('rm {0}'.format(result))
save(result, '{0} calculations'.format(id))
save(result, 'B,C,C44')

l = 3.602
atoms = bulk('Fe', 'fcc', a=l)
atoms.set_tags([1])

cr = 0.15
ni = 0.15
fe = 1.0-cr-ni

alloys = []
alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))
alloys.append(Alloy(1, 'Cr', cr / 2, 1.0))
alloys.append(Alloy(1, 'Cr', cr / 2, -1.0))
alloys.append(Alloy(1, 'Ni', ni / 2, 1.0))
alloys.append(Alloy(1, 'Ni', ni / 2, -1.0))

calc = EMTO()
calc.set(
    lat=2,
    kpts=[13, 13, 13], )

calc.set_alloys(alloys)

#v0, e0, B = EMTO_EOS(temp_dir + '/B', atoms, calc)

#l = (4.0 * v0) ** (1.0 / 3.0)

# print v0, e0, B

#C = EC_FCC_C(temp_dir + '/C', atoms, calc)

l = 3.602
atoms = bulk('Fe', 'fcc', a=l)
atoms.set_tags([1])
C44 = EC_FCC_C44(temp_dir + '/C44', atoms, calc)

save(result, '{0},{1},{2}'.format(0,C,C44))









