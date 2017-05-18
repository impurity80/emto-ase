

# *******************************************************
#  Copyright (C) Korea Institute of Materials Science - All Rights Reserved
#  Unauthorized copying of this file, via any medium is strictly prohibited
#  Proprietary and confidential
#  Written by Jae Hoon Jang <jhjang@kims.re.kr>, 18th, May, 2017
#  Version : 0.2
# *******************************************************/

import numpy as np
import os
from ase.utils.eos import EquationOfState
from ase.units import *

from emto import *

def EMTO_EOS(work_dir, atoms, calc, OPTIONS = np.linspace(0.98, 1.02, 9)):

    os.system('mkdir -p {0}/calc'.format(work_dir))

    result = '{0}/result.txt'.format(work_dir)
    os.system('rm {0}'.format(result))

    save(result, 'Bulk modululs calculation')

    volumes = []
    energies = []

    for o in OPTIONS:

        a = atoms.copy()
        scale = [[o, 0, 0], [0, o, 0], [0, 0, o]]
        a.set_cell(np.dot(atoms.get_cell(), scale), scale_atoms=True)

        calc.set(dir='{0}/calc/opt-{1:0.3f}'.format(work_dir, o))
        a.set_calculator(calc)

        n = a.get_number_of_atoms()
        e = a.get_potential_energy() / n
        v = a.get_volume() / n

        if e < -0.001:
            volumes.append(v)
            energies.append(e)

        save(result, '{0},{1},{2}'.format(o, v, e) )

    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    eos.plot('{0}/B.png'.format(work_dir))

    save(result, 'calculation finshised')

    save(result, OPTIONS)
    save(result, volumes)
    save(result, energies)

    save(result, '------------------------')

    save(result, '{0},{1},{2}', v0, e0, B/kJ*1.0e24)

    save(result, '------------------------')

    return v0, e0, B/kJ*1.0e24