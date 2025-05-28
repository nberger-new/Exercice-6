import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import trapezoid as tpz
from scipy.optimize import curve_fit
import os

plt.rc('font', family='serif')

filename = 'output_V0=1330_obs.out'  # Sostituisci con il path corretto se serve

# Carico i dati
data = np.loadtxt(filename)
fs = 17  # font size
time = data[:, 0]
prob_left = data[:, 1]
prob_right = data[:, 2]

# Plot probabilità a sinistra e destra nello stesso grafico
plt.figure(figsize=(8, 5))
plt.scatter(time, prob_left, label=r"$P(t)_{x<0}$", color='mediumvioletred', marker='o')
plt.scatter(time, prob_right, label=r"$P(t)_{x>0}$", color='mediumturquoise', marker ='o')
plt.xlabel(r"Temps [s]", fontsize=fs)
plt.ylabel(r"P(t)", fontsize=fs)
plt.grid(True)
plt.tick_params(labelsize=fs)
plt.tight_layout()
plt.show()


