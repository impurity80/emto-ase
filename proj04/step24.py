
import os
from emto import *
from emto_sfe import *

id = '24'

curr_dir = os.getcwd()
temp_dir = curr_dir.replace('work', 'temp')  + '/' + id

result = '{0}/result-{1}.txt'.format(curr_dir.replace('work', 'temp'),id)
os.system('rm {0}'.format(result))
save(result, '{0} calculations'.format(id))
# save(result, 'v0,e0,B,C,C44,C11,C12,A,Gv,Gr,Gh,C12-C44,B/Gh,Eh,Ph')

x = 0.002
mn = 0.20
fe = 1.0-x-mn

alloys = []
alloys.append(Alloy(1, 'Fe', fe / 2.0 , 1.0))
alloys.append(Alloy(1, 'Fe', fe / 2.0 , -1.0))
alloys.append(Alloy(1, 'Al', x / 2.0 , 1.0))
alloys.append(Alloy(1, 'Al', x / 2.0 , -1.0))
alloys.append(Alloy(1, 'Mn', mn / 2.0 , 1.0))
alloys.append(Alloy(1, 'Mn', mn / 2.0 , -1.0))

r = EMTO_SFE(temp_dir, alloys, 3.60)

print r

save(result, '{0}'.format(r))









