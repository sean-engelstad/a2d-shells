import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import niceplots

# r/t of each cylinder design
df = pd.read_csv("../kdfs.csv", skiprows=1)
print(df)
print(df.columns[1])
nelems = df[["nelems"]].to_numpy()
rt_values = df[[" r/t"]].to_numpy()
nasa_kdfs = df[[" nasaKDF"]].to_numpy()
tacs_kdfs = df[[" tacsKDF"]].to_numpy()

plt.style.use(niceplots.get_style())

my_colors = ["#20a39eff", "#ffba49ff", "#ef5b5bff"]

plt.vlines(x=100, ymin=0.0, ymax=1.2, color='tab:grey', linestyles='--', linewidth=2)

nelems_unique = np.unique(nelems)
mask = np.logical_and(nelems == 10000, rt_values > 10)
plt.plot(rt_values[mask], nasa_kdfs[mask], "ko--", label="nasa-kdf")
for i_nelems, c_nelems in enumerate(nelems_unique):
    nelems_mask = np.logical_and(nelems == c_nelems, rt_values > 10)
    plt.plot(rt_values[nelems_mask], tacs_kdfs[nelems_mask], "o-", color=my_colors[i_nelems], label=f"tacs-{c_nelems/1000:.0f}k")

plt.margins(x=0.05, y=0.05)
plt.xlabel("r/t")
plt.xscale('log')
plt.legend()
plt.ylabel("knockdown factor (KDF)")
plt.show()
