import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the csv file "load_disp.csv" and plot it
df = pd.read_csv("load-disp.csv")
disp = df[["|u|"]].to_numpy()
load = df[["lambda"]].to_numpy()

freq = 10
plt.plot(disp[0::freq], load[0::freq], "ko-", linewidth=2)
plt.xlabel("|u|")
plt.ylabel("lambda")
plt.show()
