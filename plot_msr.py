import matplotlib.pyplot as plt
import numpy as np
data1 = []
data2 = []
data3 = []
with open('max_sr_qpsk','r') as data:
    for line in data:
        data1.append(int(line))


data1 = np.array(data1)
smf= np.arange(30,50,0.01)
print(data1[smf==40])


with open('pscf_msr') as data:
    for line in data:
        data2.append(int(line))
print(data1[smf==40])

with open('nzdsf_msr') as data:
    for line in data:
        data3.append(int(line))
# plt.plot(np.arange(30,50,0.01),data2,label='pscf')
#
# plt.plot(np.arange(30,50,0.01),data3,label='nzdsf')
#
# plt.xlabel(r'$\Delta f_{ch} [GHZ]$')
# plt.ylabel('$N_s^{max}$')
# plt.xlim(30,50)
# plt.ylim(10,100)
# plt.legend()
# plt.show()