
import os
from ase import Atom, Atoms
from emto import *
from emto_eos import *

from ase.lattice import bulk

id = '1'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp') + '/' + id

result = '{0}/result-{1}.txt'.format(temp_dir,id)
os.system('rm {0}'.format(result))

l = 4.50
atoms = Atoms('Nb4',
              scaled_positions=[
                                (0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[l,l,l],
              pbc=(1,1,1))

atoms.set_tags([1, 1, 1, 1])

atoms = atoms + Atom('C', position=(0, 0, 0.5 * l), tag=2)


X_OPT = np.linspace(0.0, 1.0, 11)

for x in X_OPT:

    alloys = []
    alloys.append(Alloy(1, 'Nb', 1.0-x, 1.0))
    alloys.append(Alloy(1, 'Mo', x, 1.0))
    alloys.append(Alloy(2, 'C', 1.0, 0.0))

    calc = EMTO()
    calc.set(
        lat=1,
        kpts=[13, 13, 13],
        aw=0.60,
        dmax=1.4,
        amix=0.01 )

    calc.set_alloys(alloys)

    r = EMTO_EOS(temp_dir + '/B' + '/{0:0.3f}'.format(x) , atoms, calc)

    save(result, '{0}'.format(r))













