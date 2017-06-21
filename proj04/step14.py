
import os
from emto import *
from emto_sfe import *

id = '14'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)
os.system('rm {0}'.format(result))
save(result, '{0} calculations'.format(id))
# save(result, 'v0,e0,B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44,B/Gh,Eh,Ph')

cr = 0.15
fe = 1.0-cr

alloys = []
alloys.append(Alloy(1, 'Fe', fe / 2.0 , 1.0))
alloys.append(Alloy(1, 'Fe', fe / 2.0 , -1.0))
alloys.append(Alloy(1, 'Cr', cr / 2.0 , 1.0))
alloys.append(Alloy(1, 'Cr', cr / 2.0 , -1.0))

r = EMTO_SFE(temp_dir, alloys, 3.60)

print r

save(result, '{0}'.format(r))









