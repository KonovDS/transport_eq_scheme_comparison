import matplotlib.pyplot as plt
import numpy as np

f0, l0 = open("./0.txt", "r"), []
f1, l1 = open("./1.txt", "r"), []
f2, l2 = open("./2.txt", "r"), []
f3, l3 = open("./3.txt", "r"), []
f4, l4 = open("./4.txt", "r"), []
f5, l5 = open("./5.txt", "r"), []

for line in f0:   
    for w in line.split():         
        l0.append(float(w))

for line in f1:   
    for w in line.split():         
        l1.append(float(w))

for line in f2:   
    for w in line.split():         
        l2.append(float(w))

for line in f3:   
    for w in line.split():         
        l3.append(float(w))

for line in f4:   
    for w in line.split():         
        l4.append(float(w))

for line in f5:   
    for w in line.split():         
        l5.append(float(w))

fig, ax = plt.subplots()

ax.plot(l0, label='Initial condition')
ax.plot(l1, label='Godunov scheme')
ax.plot(l2, label='MacCormack scheme')
ax.plot(l3, label='Holodnov scheme')
ax.plot(l4, label='Cubic approximation scheme')
ax.plot(l5, label='Hybrid scheme')
ax.legend()

fig.set_figheight(5)
fig.set_figwidth(8)
plt.show()