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
import matplotlib.pyplot as plt

# Example data
x = [1, 2, 3, 4, 5]
y1 = [1, 4, 9, 16, 25]
y2 = [1, 2, 3, 4, 5]
y3 = [25, 20, 15, 10, 5]

# Plot multiple series
plt.plot(x, y1, label='y = x^2')
plt.plot(x, y2, label='y = x')
plt.plot(x, y3, label='y = 30 - 5x')

# Add title and labels
plt.title('Multiple Series on the Same Plot')
plt.xlabel('x-axis')
plt.ylabel('y-axis')

# Add a legend
plt.legend()

# Show the plot
plt.show()
