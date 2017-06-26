
import os
from emto import *
from emto_sfe import *

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

id = '.'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)

if rank==0:
    os.system('rm {0}'.format(result))
    save(result, '{0} calculations'.format(id))

X_OPT = np.linspace(0.0, 0.15, 16)

x = X_OPT[rank]
mn = 0.20
fe = 1-x-mn

alloys = []
alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))
alloys.append(Alloy(1, 'Si', x / 2, 1.0))
alloys.append(Alloy(1, 'Si', x / 2, -1.0))
alloys.append(Alloy(1, 'Mn', mn / 2.0 , 1.0))
alloys.append(Alloy(1, 'Mn', mn / 2.0 , -1.0))

r = EMTO_SFE(temp_dir+'/'+'{0:0.2f}'.format(x), alloys, 3.60)

print r

save(result, '{0}, {1}'.format(x, r))
