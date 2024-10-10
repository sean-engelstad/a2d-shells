import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8, 6))  # 2 rows, 1 column

w = np.linspace(0.0, 2.5, 100)
ax2.plot(w, 0.0 * w, "k--")

for i, KS in enumerate([0.0, 3.0, 5.0]):
    # KS = 0.0 # 0.0, 3.0, 5.0 KS values
    E = 10e7; A = 0.1
    L = np.sqrt(100**2 + 1)
    k = E * A * L
    strain = (np.sqrt(100**2 + (1 - w)**2) - L) / L
    U = 0.5 * k * strain**2
    U2 = 0.5 * KS * w**2
    # Pi = U - P * w
    # equilibrium when dPi/dw = 0 => solve for P
    dstrain_dw = (w-1) / np.sqrt(100**2 + (1-w)**2) / L
    dU2dw = KS * w
    P = k * strain * dstrain_dw + dU2dw
    ax1.plot(w, P, label=f"{KS=}")
    ax1.set_xlabel(f"w disp")
    ax1.set_ylabel(f"P load")

    # also compute the 2nd derivatives for stability
    dstrain_dww = 1.0 / np.sqrt(100**2 + (1-w)**2) / L - 1.0 * dstrain_dw / np.sqrt(100**2 + (1-w)**2)
    stability = k * dstrain_dw**2 + k * strain * dstrain_dww + KS
    # print(f"{stability=}")
    ax2.plot(w, stability, "--", label=f"stab{i}")
    ax2.set_xlabel("w disp")
    ax2.set_ylabel("stability d^2Pi/dw^2")

plt.legend()
plt.show()