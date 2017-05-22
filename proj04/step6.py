
import os
from emto import *
from emto_eos import *
from ec_fcc import *

from ase.lattice import bulk

id = '6'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)
os.system('rm {0}'.format(result))
save(result, '{0} calculations'.format(id))
save(result, 'v0,e0,B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44')

ni = 0.074
al = 0.30
mn = 0.34

fe = 1.0-ni-mn-al

alloys = []
alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))
alloys.append(Alloy(1, 'Mn', mn / 2, 1.0))
alloys.append(Alloy(1, 'Mn', mn / 2, -1.0))
alloys.append(Alloy(1, 'Al', al / 2, 1.0))
alloys.append(Alloy(1, 'Al', al / 2, -1.0))
alloys.append(Alloy(1, 'Ni', ni / 2, 1.0))
alloys.append(Alloy(1, 'Ni', ni / 2, -1.0))

v0, e0, B, C, C44, C11, C12, A, Gv, Gr, Gh = fcc_elastic_constants(temp_dir, alloys, 3.7)

save(result, '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(v0, e0, B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44))









