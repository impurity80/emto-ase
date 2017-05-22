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
import matplotlib.pyplot as plt
from ase.lattice import bulk

from emto_eos import *
from emto import *

import numpy.polynomial.polynomial as poly

def bcc_elastic_constants(work_dir, alloys, l=2.86):
    atoms = bulk('Fe', 'bcc', a=l)
    atoms.set_tags([1])

    calc = EMTO()
    calc.set(
        lat=3,
        kpts=[13, 13, 13],
        dmax=2.20,)

    calc.set_alloys(alloys)

    v0, e0, B = EMTO_EOS(work_dir + '/B', atoms, calc)

    print v0, e0, B

    l = (2.0 * v0) ** (1.0 / 3.0)
    atoms = bulk('Fe', 'bcc', a=l)
    atoms.set_tags([1])

    C = EC_BCC_C(work_dir + '/C', atoms, calc)
    C44 = EC_BCC_C44(work_dir + '/C44', atoms, calc)
    C11 = (3.0 * B + 2.0 * C) / 3.0
    C12 = (3.0 * B - C) / 3.0
    A = 2.0 * C44 / C
    Gv = (C + 3.0 * C44) / 5.0
    Gr = (5.0 * C * C44) / (4.0 * C44 + 3.0 * C)
    Gh = (Gv + Gr) / 2.0

    return v0, e0, B, C, C44, C11, C12, A, Gv, Gr, Gh


def EC_BCC_C(work_dir, atoms, calc, OPTIONS = np.linspace(0.00, 0.05, 6)):

    os.system('mkdir -p {0}/calc'.format(work_dir))

    result = '{0}/result.txt'.format(work_dir)
    os.system('rm {0}'.format(result))

    save(result, 'FCC C11-C12 calculation')

    volumes = []
    energies = []

    for o in OPTIONS:

        a = atoms.copy()
        scale = [[1 + o, 0, 0], [0, 1 - o, 0], [0, 0, 1 / (1 - o ** 2)]]
        a.set_cell(np.dot(atoms.get_cell(), scale), scale_atoms=True)

        calc.set(dir='{0}/calc/opt-{1:0.3f}'.format(work_dir, o),
                 lat=10,
                 kpts=[27, 27, 27],
                 dmax=2.20,)

        a.set_calculator(calc)

        n = a.get_number_of_atoms()
        e = a.get_potential_energy() / n
        v = a.get_volume() / n

        if e < -0.001:
            volumes.append(v)
            energies.append(e)

        save(result, '{0},{1},{2}'.format(o, v, e) )

    OPT = (-1.0 * OPTIONS[::-1]).tolist() + OPTIONS.tolist()
    energies = (energies[::-1]) + energies

    coefs = poly.polyfit(OPT, energies, 3)

    C = coefs[2] / volumes[0] / kJ * 1.0e24

    print C

    x_new = np.linspace(OPT[0] - 0.01, OPT[-1] + 0.01, num=len(OPT) * 10)

    ffit = poly.polyval(x_new, coefs)

    plt.scatter(OPT, energies)
    plt.plot(x_new, ffit)
    plt.savefig('C.png')

    os.system('mv C.png {0}'.format(work_dir))

    save(result, OPT)
    save(result, volumes)
    save(result, energies)

    save(result, '------------------------')

    save(result, '{0}'.format(C))

    save(result, '------------------------')

    return C

def EC_BCC_C44(work_dir, atoms, calc, OPTIONS = np.linspace(0.00, 0.05, 6)):

    os.system('mkdir -p {0}/calc'.format(work_dir))

    result = '{0}/result.txt'.format(work_dir)
    os.system('rm {0}'.format(result))

    save(result, 'FCC C44 calculation')

    volumes = []
    energies = []

    for o in OPTIONS:

        cell = atoms.get_cell()
        l = np.linalg.norm(cell[0])

        a = bulk('Fe', 'fcc', a=l)
        a.set_tags(atoms.get_tags())

        dist = [[np.sqrt(2)*(1 + o), 0, 0], [0, np.sqrt(2)*(1 - o), 0], [0, 0, 1.0 / (1 - o ** 2)]]
        a.set_cell(np.dot(a.get_cell(), dist), scale_atoms=True)

        calc.set(dir='{0}/calc/opt-{1:0.3f}'.format(work_dir, o),
                 lat=11,
                 dmax=1.60)

        a.set_calculator(calc)

        n = a.get_number_of_atoms()
        e = a.get_potential_energy() / n
        v = a.get_volume() / n

        if e < -0.001:
            volumes.append(v)
            energies.append(e)

        save(result, '{0},{1},{2}'.format(o, v, e) )

    OPT = (-1.0 * OPTIONS[::-1]).tolist() + OPTIONS.tolist()
    energies = (energies[::-1]) + energies

    coefs = poly.polyfit(OPT, energies, 3)

    C44 = coefs[2] / volumes[0] / kJ * 1.0e24 / 2.0

    x_new = np.linspace(OPT[0] - 0.01, OPT[-1] + 0.01, num=len(OPT) * 10)

    ffit = poly.polyval(x_new, coefs)

    plt.scatter(OPT, energies)
    plt.plot(x_new, ffit)
    plt.savefig('C44.png')

    os.system('mv C.png {0}'.format(work_dir))

    save(result, OPT)
    save(result, volumes)
    save(result, energies)

    save(result, '------------------------')

    save(result, '{0}'.format(C44))

    save(result, '------------------------')

    return C44