import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False

file_list = [
    "output_type2_V0=-3000_obs.out"
    
]

potentials_list = [
    0, 10, 20, 40, 100, 300
]

V0_list = []
E_over_V0 = []
Ptrans_list = []

for filename, V0 in zip(file_list, potentials_list):
    try:
        data = np.loadtxt(filename)
        t_index = np.argmin(np.abs(data[:, 0] - 0.035))
        P_right = data[t_index, 2]
        E = data[t_index, 3]

        V0_list.append(V0)
        E_over_V0.append(E / V0)
        Ptrans_list.append(P_right)
    except Exception as e:
        print(f"Errore nel file {filename}: {e}")

E_over_V0 = np.array(E_over_V0)
Ptrans_list = np.array(Ptrans_list)

# Ordina i dati per E_over_V0 per tracciare linee ordinate
sorted_indices = np.argsort(E_over_V0)
E_over_V0 = E_over_V0[sorted_indices]
Ptrans_list = Ptrans_list[sorted_indices]
V0_list = np.array(V0_list)[sorted_indices]

plt.figure(figsize=(10, 6))

# Aggiungi zone colorate orizzontali per le probabilità
plt.axhspan(0.7, 1.0, color='lightgreen', alpha=0.3, label='Probabilità > 0.7')
plt.axhspan(0.5, 0.7, color='khaki', alpha=0.3, label='Probabilità tra 0.5 e 0.7')
plt.axhspan(0.0, 0.5, color='lightcoral', alpha=0.3, label='Probabilità < 0.5')

# Traccia linee per ogni potenziale, usando colori diversi (qui semplicemente cycling)
colors = plt.cm.viridis(np.linspace(0, 1, len(V0_list)))

for i, (V0, color) in enumerate(zip(V0_list, colors)):
    plt.plot(E_over_V0[i:i+2], Ptrans_list[i:i+2], '-', color=color, label=f"$V_0={V0}$")  # Questo però traccia segmenti (a meno che non ci siano 2 punti per potenziale)

# Poiché per ogni potenziale hai solo un punto (E_over_V0, Ptrans), ha senso tracciare i punti collegati
# Se invece vuoi tracciare solo i punti con linee, si potrebbe fare così:

plt.plot(E_over_V0, Ptrans_list, '-o', color='coral', label=r"$P_{\mathrm{trans}}$")

plt.xlabel(r"$\langle E \rangle / V_0$", fontsize=20)
plt.ylabel(r"P$_{\mathrm{trans}}$", fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.grid(True)

# Mostra legenda per zone e per curva
plt.tight_layout()
plt.show()
