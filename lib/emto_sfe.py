
# *******************************************************
#  Copyright (C) Korea Institute of Materials Science - All Rights Reserved
#  Unauthorized copying of this file, via any medium is strictly prohibited
#  Proprietary and confidential
#  Written by Jae Hoon Jang <jhjang@kims.re.kr>, 18th, May, 2017
#  Version : 0.3
# *******************************************************/

import numpy as np
import os
from ase.utils.eos import EquationOfState
from ase.units import *
from ase.lattice import bulk
from ase import Atom, Atoms

from emto_eos import *
from emto import *

def EMTO_SFE(work_dir, alloys, l=3.60):

    atoms = bulk('Fe', 'fcc', a=l)
    atoms.set_tags([1])

    calc = EMTO()
    calc.set(
        lat=2,
        kpts=[17, 17, 17],
        amix=0.02, )

    calc.set_alloys(alloys)

    v0, e0, B0 = EMTO_EOS(work_dir + '/FCC', atoms, calc, np.linspace(0.98, 1.02, 9))
    l = (4.0 * v0) ** (1.0 / 3.0)
    a0 = l / np.sqrt(2)
    c0 = np.sqrt(8 / 3.0) * a0

    os.system('mkdir -p {0}/SFE'.format(work_dir))

    result = '{0}/SFE/result.txt'.format(work_dir)
    os.system('rm {0}'.format(result))

    save(result, 'sfe calculation')

    # HCP ab sequence
    atoms = Atoms('Fe2',
                  scaled_positions=[(0, 0, 0),
                                    (1. / 3., 1. / 3., 1. / 2.)],
                  cell=[[1. / 2., sqrt(3) / 2., 0], [-1. / 2., sqrt(3) / 2., 0], [0, 0, 1.0]],
                  pbc=(1, 1, 1))

    scale = [[a0, 0, 0], [0, a0, 0], [0, 0, c0]]
    atoms.set_cell(np.dot(atoms.get_cell(), scale), scale_atoms=True)

    atoms.set_tags([1, 1])

    calc.set(dir='{0}/SFE/ab'.format(work_dir),
             lat=9,
             kpts=[13,13,11],
             dmax=2.4,
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)

    e1 = atoms.get_potential_energy() / atoms.get_number_of_atoms()
    v1 = atoms.get_volume() / atoms.get_number_of_atoms()

    save(result, 'ab sequence : {0},{1},{2}'.format(l, v1, e1))

    # FCC abc sequence
    atoms = Atoms('Fe3',
                  scaled_positions=[(0, 0, 0),
                                    (1. / 3., 1. / 3., 1. / 3.),
                                    (2. / 3., 2. / 3., 2. / 3.)],
                  cell=[[1. / 2., sqrt(3) / 2., 0], [-1. / 2., sqrt(3) / 2., 0], [0, 0, 1.5]],
                  pbc=(1, 1, 1))

    scale = [[a0, 0, 0], [0, a0, 0], [0, 0, c0]]
    atoms.set_cell(np.dot(atoms.get_cell(), scale), scale_atoms=True)

    atoms.set_tags([1, 1, 1])

    calc.set(dir='{0}/SFE/abc'.format(work_dir),
             lat=9,
             kpts=[13,13,9],
             dmax=2.4,
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)

    e2 = atoms.get_potential_energy() / atoms.get_number_of_atoms()
    v2 = atoms.get_volume() / atoms.get_number_of_atoms()

    save(result, 'abc sequence : {0},{1},{2}'.format(l, v2, e2))

    # double hcp abac sequence
    atoms = Atoms('Fe4',
                  scaled_positions=[(0, 0, 0),
                                    (1. / 3., 1. / 3., 1. / 4.),
                                    (0, 0, 2./4.),
                                    (2. / 3., 2. / 3., 3. / 4.)],
                  cell=[[1. / 2., sqrt(3) / 2., 0], [-1. / 2., sqrt(3) / 2., 0], [0, 0, 2.0]],
                  pbc=(1, 1, 1))

    scale = [[a0, 0, 0], [0, a0, 0], [0, 0, c0]]
    atoms.set_cell(np.dot(atoms.get_cell(), scale), scale_atoms=True)

    atoms.set_tags([1, 1, 1, 1])

    calc.set(dir='{0}/SFE/abac'.format(work_dir),
             lat=9,
             kpts=[13,13,7],
             dmax=2.4,
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)

    e3 = atoms.get_potential_energy() / atoms.get_number_of_atoms()
    v3 = atoms.get_volume() / atoms.get_number_of_atoms()

    save(result, 'abac sequence : {0},{1},{2}'.format(l, v3, e3))

    area = sqrt(3)/2.0*a0*a0
    unit = 1.60210*10*1000 # eV/A^2 -> mJ/m^2

    return l, a0, c0, B0, v0, e0, e1, e2, e3, area , 2.0*(e1-e2), e1+2.0*e3-3.0*e2, 2.0*(e1-e2)/area*unit,  (e1+2.0*e3-3.0*e2)/area*unit