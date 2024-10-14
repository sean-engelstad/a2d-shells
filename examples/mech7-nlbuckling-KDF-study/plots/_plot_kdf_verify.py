import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# r/t of each cylinder design
rt_values = np.array([1000, 500, 300, 100, 50, 25, 10])
nasa_kdfs = np.array([0.223, 0.322, 0.404, 0.581, 0.678, 0.758, 0.82])
# from 10k element data
tacs_kdfs = np.array([0.61, 0.664, 0.677, 0.691, 0.705, 0.763, 0.82])

plt.plot(rt_values[::-1], nasa_kdfs[::-1], "o-", label="nasa-kdf")
plt.plot(rt_values[::-1], tacs_kdfs[::-1], "o-", label="tacs-mode1-kdf")

plt.xlabel("r/t")
plt.xscale('log')
plt.legend()
plt.ylabel("knockdown factor (KDF)")
plt.show()
