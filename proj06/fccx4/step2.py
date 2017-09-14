
import os
from ase import Atom, Atoms
from emto import *
from emto_eos import *

from ase.lattice import bulk

id = '2'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp') + '/' + id

result = '{0}/result-{1}.txt'.format(temp_dir,id)
os.system('rm {0}'.format(result))

l = 3.70
atoms = Atoms('Fe4',
              scaled_positions=[
                                (0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[l,l,l],
              pbc=(1,1,1))

atoms.set_tags([1, 1, 1, 1])

atoms = atoms + Atom('C', position=(0, 0, 0.5 * l), tag=2)

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
alloys.append(Alloy(2, 'C', 0.5, 0.0))
alloys.append(Alloy(2, 'Va', 0.5, 0.0))

calc = EMTO()
calc.set(
    lat=1,
    kpts=[13, 13, 13],
    aw=0.60,
    dmax=1.0,
    amix=0.02 )

calc.set_alloys(alloys)

r = EMTO_EOS(temp_dir + '/B', atoms, calc)

save(result, '{0}'.format(r))













