import numpy as np
import matplotlib.pyplot as plt

# Impostazioni grafico
plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False

# Lista file e tipi di potenziale
file_list = ["output1_obs.out", "output2_obs.out", "output3_obs.out", "output4_obs.out"]
potentials_types = [1, 2, 3, 4]
colors = ['royalblue', 'salmon', 'gold', 'springgreen']

# Creazione figura per il grafico
plt.figure(figsize=(8, 5))

# Ciclo su ogni file

for filename, potential, color in zip(file_list, potentials_types, colors):
    try:
        data = np.loadtxt(filename)
        t = data[:, 0]
        x_moy = data[:, 4]
        x2_moy = data[:, 5]
        delta_x = np.sqrt(x2_moy - x_moy**2)

        # Colore personalizzato per ogni curva
        plt.plot(t, delta_x, label=f'Potentiel {potential}', color=color)

        print(f"[Potentiel {potential}] Moyenne de Δx : {np.mean(delta_x):.5f}")
        print(f"[Potentiel {potential}] Écart-type de Δx : {np.std(delta_x):.5e}")

    except Exception as e:
        print(f"Erreur dans le fichier {filename} : {e}")
        continue

# Etichette e layout
plt.xlabel("Temps [s]", fontsize = 20)
plt.ylabel(r"$\Delta x$", fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.grid(True)
plt.legend(fontsize=17)
plt.tight_layout()
plt.show()
