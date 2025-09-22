#numpytest.py
import numpy as np


m=np.zeros((3,4))
print(m)

temps = list(range(300,2001,100))
temps.insert(0,298.15)
print(temps)

for i in range(6):
    print(i)

s = ""
for i in range(3):
    for j in range(4):
        m[i][j]=str((i+1)*(j+1))
        s += str((i+1)*(j+1))
        s+=","
    
    s+="\n"

print(s)

print(m)

f = ""

for j in range(4):
    for i in range(3):
        f+=str("{:.2f}".format(m[i][j]))
        f+=","
    f+="\n"

print(f)