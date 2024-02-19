import numpy as np 
import matplotlib.pyplot as plt

###############################
#  POSTPROCESSING & PLOTTING  #
###############################

x = np.genfromtxt("dump.txt")

x0 = x[:, 0]
x1 = x[:, 1]

kern_half_size = 5

kern = list(range(1, kern_half_size+2))
kern.extend(list(range(kern_half_size, 0, -1)))

kern = np.array(kern)
kern = kern / np.sum(kern)

def extend(x):
    return np.hstack((
        x[-kern_half_size:],
        x,
        x[:kern_half_size]
        ))

x0_extended = extend(x0)
x1_extended = extend(x1)

y0 = np.convolve(x0_extended, kern, mode="valid")
y1 = np.convolve(x1_extended, kern, mode="valid")

with open("out.csv", "w") as out:
    for i, y in enumerate(y0):
        out.write("%f,%f\n" % (y, y1[i]))

print("Result writting to out.csv")

plt.plot(x0, x1)
plt.plot(y0, y1)

plt.show()
