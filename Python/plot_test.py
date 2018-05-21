import numpy as np
import matplotlib.pyplot as plt


a = np.array([1,2,3,6,7])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(a,label=str(0.01))
ax.plot(a-2,label="goodbye")
ax.legend()
plt.show()