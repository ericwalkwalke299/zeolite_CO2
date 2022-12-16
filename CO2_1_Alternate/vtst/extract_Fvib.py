import os
import re
import numpy as np

T = input("Enter Temp (K): ") 
T = float(T)
os.system('/projects/academic/mdupuis2/software/vtst/vtsttools/vtstscripts/dymmatrix.pl')
freq = []
with open('freq.dat') as f:
    content = f.readlines()
    for line in content:
        line = line.strip()
        nums = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        freq.append(float(nums[0]))
freq = np.array(freq)
l = len(freq)
for k in range(l):
    if freq[k] < 100:
        freq[k] = 100
c = 30000000000
kb = 0.00008617
h = 0.00000000000000413566
freq = freq * c
zpe = .5*h*freq
ZPE = np.sum(zpe)
qvib = 1/(1-np.exp((-h*freq)/(kb*T)))
Qvib = np.prod(qvib)
Fvib = ZPE + kb*T*np.log(1/Qvib)
S_T = kb*T*np.log(Qvib)
print(ZPE)
print(Qvib)
print(Fvib)
print(S_T)
