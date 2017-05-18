
import os
from emto import *
from emto_eos import *
from ec_fcc import *

from ase.lattice import bulk

id = '4'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)
os.system('rm {0}'.format(result))
save(result, '{0} calculations'.format(id))
save(result, 'B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44')

l = 3.602
atoms = bulk('Fe', 'fcc', a=l)
atoms.set_tags([1])

ni = 0.28
co = 0.17
al = 0.114

fe = 1.0-ni-co-al

alloys = []
alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))
alloys.append(Alloy(1, 'Co', co / 2, 1.0))
alloys.append(Alloy(1, 'Co', co / 2, -1.0))
alloys.append(Alloy(1, 'Al', al / 2, 1.0))
alloys.append(Alloy(1, 'Al', al / 2, -1.0))
alloys.append(Alloy(1, 'Ni', ni / 2, 1.0))
alloys.append(Alloy(1, 'Ni', ni / 2, -1.0))

calc = EMTO()
calc.set(
    lat=2,
    kpts=[13, 13, 13], )

calc.set_alloys(alloys)

v0, e0, B = EMTO_EOS(temp_dir + '/B', atoms, calc)

l = (4.0 * v0) ** (1.0 / 3.0)
atoms = bulk('Fe', 'fcc', a=l)
atoms.set_tags([1])

C = EC_FCC_C(temp_dir + '/C', atoms, calc)
C44 = EC_FCC_C44(temp_dir + '/C44', atoms, calc)
C11 = (3.0*B+2.0*C)/3.0
C12 = (3.0*B-C)/3.0
A = 2.0*C44/C
Gv = (C+3.0*C44)/5.0
Gr = (5.0*C*C44)/(4.0*C44+3.0*C)
Gh = (Gv+Gr)/2.0

save(result, '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}'.format(B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44))

