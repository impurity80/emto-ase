
import os
from emto import *
from ec_fcc import *

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

id = '8'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)

if rank==0:
    os.system('rm {0}'.format(result))
    save(result, '{0} calculations'.format(id))
    save(result, 'ni, cr, v0,e0,B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44')

X_OPT = np.linspace(0.08, 0.23, 16)
Y_OPT = np.linspace(0.13, 0.23, 11)

ni = X_OPT[rank]
for cr in Y_OPT:
    fe = 1.0-cr-ni

    alloys = []
    alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))
    alloys.append(Alloy(1, 'Cr', cr / 2, 1.0))
    alloys.append(Alloy(1, 'Cr', cr / 2, -1.0))
    alloys.append(Alloy(1, 'Ni', ni / 2, 1.0))
    alloys.append(Alloy(1, 'Ni', ni / 2, -1.0))

    v0, e0, B, C, C44, C11, C12, A, Gv, Gr, Gh = fcc_elastic_constants(temp_dir+'/'+'{0:0.2f}-{1:0.2f}'.format(ni, cr), alloys)

    save(result, '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}'.format(ni, cr, v0, e0, B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44))
