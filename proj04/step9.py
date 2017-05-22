
import os
from emto import *
from ec_bcc import *

id = '9'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)
os.system('rm {0}'.format(result))
save(result, '{0} calculations'.format(id))
save(result, 'v0,e0,B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44')

cr = 0.02
ni = 0.02
fe = 1.0-cr-ni

alloys = []
alloys.append(Alloy(1, 'Fe', fe , 1.0))
alloys.append(Alloy(1, 'Cr', cr , 1.0))
alloys.append(Alloy(1, 'Ni', ni , 1.0))

v0, e0, B, C, C44, C11, C12, A, Gv, Gr, Gh = bcc_elastic_constants(temp_dir, alloys)

save(result, '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}'.format(v0, e0, B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44))









