import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import trapezoid as tpz
from scipy.optimize import curve_fit
import pdb
import os
import time

plt.rc('font',family='serif')

from matplotlib import cm
from matplotlib.ticker import LinearLocator

# Variation de l'énergie : variation par rapport à la position précédente : (E[i] - E[i+1])/dt

# Parameters
executable = './Exercice6_2025_student' # Nome dell'eseguibile
repertoire = r"/home/nbh/python/Physique-Numerique/Exercice-6"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Nome del file di input

one = True
two = False
three = False
init = False

fs = 15

if init :
  Nsteps = 800
  Nintervals = 512
  V0 = [0, 1000, 2000, 3000]
  xa = [0, -0.5, -0.5, -0.5]
  xb = [0, 0.5, 0.5, 0.5]

if one or two :
  Nsteps = [800]
  Nintervals = [512]

elif three :
  time = True
  if time :
    Nintervals = [512]
    Nsteps = [100, 200, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000]
  else :
    Nintervals = [10, 50, 100, 150, 300, 500, 1000, 1500, 3000, 5000, 10000]
    Nsteps = [800]

# Simulations
outputs = {}  # Dictionary to store output file names

if (one or two or three) :
  for i in range (len(Nsteps)):
    for j in range (len(Nintervals)) :
      output_base = f"Nsteps={Nsteps[i]}_Nintervals={Nintervals[j]}" # naming the output file
      outputs[(i, j)] = output_base
      print(f"Running simulation with Nsteps={Nsteps[i]}, Nintervals={Nintervals[j]}")
      cmd = f"{executable} {input_filename} Nsteps={Nsteps[i]} Nintervals={Nintervals[j]} output={output_base}"
      print(cmd)
      subprocess.run(cmd, shell=True)
      print('Done.')

if init :
  for i in range (len(V0)) :
    output_base = f"V0={V0[i]}" # naming the output file
    outputs[i] = output_base  # Storing the input files by fixed om for future use
    print(f"Running simulation with V0 = {V0[i]}")
    cmd = f"{executable} {input_filename} V0={V0[i]} xa={xa[i]} xb={xb[i]} Nsteps={Nsteps} Nintervals={Nintervals} output={output_base}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

if init :
  fig, ax0 = plt.subplots(constrained_layout=True)
  for i, output_base in outputs.items() :
    try :
      pot_file = f'{output_base}_pot.out'
      data_pot = np.loadtxt(pot_file)

      V = data_pot[:, 1]
      xs = data_pot[:, 0]

    except :
      print('Could not open the desired files : make sure they exist !')

    color = ['royalblue', 'fuchsia', 'darkorange', 'limegreen']
    ax0.plot(xs, V, marker = 'o', color=color[i], mfc='white', label=rf'V0 = {V0[i]}, xa = {xa[i]}, xb = {xb[i]}')
    ax0.set_xlabel(r'x [m]', fontsize=fs)
    ax0.set_ylabel(r"Potential V(x) [J]", fontsize=fs)
    ax0.tick_params(axis="both", labelsize=fs)
    ax0.set_xlim(xmin = -1, xmax = 1)
    ax0.legend(fontsize=fs)
    ax0.grid(True)

else :
  for (i, j), output_base in outputs.items() :
    print(output_base)
    obs_file = f'{output_base}_obs.out'
    pot_file = f'{output_base}_pot.out'
    psi_file = f'{output_base}_psi.out'
    inc_file = f'{output_base}_inc.out'

    try :
      data_obs = np.loadtxt(obs_file)
      data_pot = np.loadtxt(pot_file)
      data_psi = np.loadtxt(psi_file)
      data_inc = np.loadtxt(inc_file)

      psi_abs = data_psi[:, 0::3]
      psi_real = data_psi[:, 1::3]
      psi_im = data_psi[:, 2::3]

      V = data_pot[:, 1]
      xs = data_pot[:, 0]

      ts = data_obs[:, 0]
      p_infs = data_obs[:, 1]
      p_sups = data_obs[:, 2]
      E = data_obs[:, 3]
      x_mean = data_obs[:, 4]
      x2_mean = data_obs[:, 5]
      p_mean = data_obs[:, 6]
      p2_mean = data_obs[:, 7]

      delta_xs = data_inc[:, 1]
      delta_ps = data_inc[:, 2]

    except :
      print('Could not open the desired files : make sure they exist !')

    hbar = 1.0546e-34

    if one :
      initial = False
      surface = False
      comparison = True

      if initial :
        fig, ax1 = plt.subplots(constrained_layout=True)
        ax1.plot(xs, psi_abs[0, :], marker = 'd', color='royalblue', mfc='white', label=rf'$|\Psi(0, x)|$, n$_x$ = {Nintervals[j]}')
        ax1.plot(xs, psi_real[0, :], marker = 'd', color='limegreen', mfc='white', label=rf'Re($\Psi(0, x)$), n$_x$ = {Nintervals[j]}')
        ax1.plot(xs, psi_im[0, :], marker = 'd', color='deeppink', mfc='white', label=rf'Im($\Psi(0, x)$), n$_x$ = {Nintervals[j]}')
        ax1.set_xlabel(r'x [m]', fontsize=fs)
        ax1.set_ylabel(r"Fonction d'onde [m]", fontsize=fs)
        ax1.tick_params(axis="both", labelsize=fs)
        ax1.legend(fontsize=fs)
        ax1.grid(True)

      if surface :
      #Surface plot as a function of (x, t)

        TT, XX = np.meshgrid(ts, xs, indexing='ij')  # TT, XX shape: (nt, nx)

        fig, ax2 = plt.subplots(constrained_layout=True)
        pcm = ax2.pcolormesh(XX, TT, psi_abs, shading='auto', cmap='plasma')
        cbar = plt.colorbar(pcm, ax=ax2)
        cbar.set_label(r"$|\Psi\,$(x, t)| [u]", fontsize=fs)
        ax2.set_xlabel("Position x [u]", fontsize=fs)
        ax2.set_ylabel(r"t/v$_0$ [u]", fontsize=fs)

        fig, ax3 = plt.subplots(constrained_layout=True)
        pcm = ax3.pcolormesh(XX, TT, psi_real, shading='auto', cmap='PuOr')
        cbar = plt.colorbar(pcm, ax=ax3)
        cbar.set_label(r"Re($\Psi\,$(x, t)) [u]", fontsize=fs)
        ax3.set_xlabel("Position x [u]", fontsize=fs)
        ax3.set_ylabel(r"t/v$_0$ [u]", fontsize=fs)

      if comparison :
        def x_th(ts_, m, x0_, omega0, p0_) :
          th = [x0_*np.cos(omega0*t) + p0_*np.sin(omega0*t)/(m*omega0) for t in ts_]
          return th

        def p_th(ts_, m, x0_, omega0, p0_) :
          th = [-m*x0_*omega0*np.sin(omega0*t) + p0_*np.cos(omega0*t) for t in ts_]
          return th

        x0 = -0.5
        om0 = 100
        p0 = p_mean[0]
        m = 1

        p_ths = p_th(ts, m, x0, om0, p0)
        x_ths = x_th(ts, m, x0, om0, p0)

        fig, ax4 = plt.subplots(constrained_layout=True)
        ax4.plot(ts, x_mean, marker = 'o', color='limegreen', mfc='white', label=r'$\langle x \rangle (t)$')
        ax4.plot(ts, x_ths, marker = 'o', color='deeppink', mfc='white', label=r'$x_{class}$(t)')
        ax4.set_xlabel(r't [s]', fontsize=fs)
        ax4.set_ylabel(r"Position [u]", fontsize=fs)
        ax4.tick_params(axis="both", labelsize=fs)
        ax4.legend(fontsize=fs)
        ax4.grid(True)

        fig, ax5 = plt.subplots(constrained_layout=True)
        ax5.plot(ts, p_mean, marker = 'o', color='limegreen', mfc='white', label=r'$\langle p \rangle (t)$')
        ax5.plot(ts, p_ths, marker = 'o', color='deeppink', mfc='white', label=r'$p_{class}$(t)')
        ax5.set_xlabel(r't [s]', fontsize=fs)
        ax5.set_ylabel(r"Quantité de mouvement [kg$\cdot$u$\cdot$s$^{-1}$]", fontsize=fs)
        ax5.tick_params(axis="both", labelsize=fs)
        ax5.legend(fontsize=fs)
        ax5.grid(True)

    if two :
      #Probabilite totale reste toujours egale a 1
      proba = False
      #Enegie reste constante
      energy = True
      uncert = True
      if proba :
        p_tot = [p_inf + p_sup for p_inf, p_sup in zip(p_infs, p_sups)]
        fig, ax9 = plt.subplots(constrained_layout=True)
        ax9.plot(ts, p_tot, marker = 'o', color='darkorange', mfc='white', label=rf'N_{{steps}} = {Nsteps[j]}')
        ax9.set_xlabel(r't [s]', fontsize=fs)
        ax9.set_ylabel(r"Probabilité totale $P_{tot}$", fontsize=fs)
        ax9.tick_params(axis="both", labelsize=fs)
        ax9.legend(fontsize=fs)
        ax9.grid(True)

      if energy :
        fig, ax6 = plt.subplots(constrained_layout=True)
        ax6.plot(ts, E, marker = 'o', color='royalblue', mfc='white', label=rf'N_{{steps}} = {Nsteps[j]}')
        ax6.set_xlabel(r't [s]', fontsize=fs)
        ax6.set_ylabel(r"Hamiltonien [J]", fontsize=fs)
        ax6.tick_params(axis="both", labelsize=fs)
        ax6.legend(fontsize=fs)
        ax6.grid(True)

      if uncert:
        # inc_xs = [np.sqrt(x2_mean_ - x_mean_**2) for x2_mean_, x_mean_ in zip(x2_mean, x_mean)]
        # inc_ps = [np.sqrt(p2_mean_ - p_mean_**2) for p2_mean_, p_mean_ in zip(p2_mean, p_mean)]
        # res_x = delta_xs == inc_xs
        # res_p = delta_ps == inc_ps
        # print(f'Res x = {res_x}, Res p = {res_p}')
        # inc_prod = [delta_x*delta_p for delta_x, delta_p in zip(inc_xs, inc_ps)]

        delta_prod = [delta_x*delta_p for delta_x, delta_p in zip(delta_xs, delta_ps)]

        fig, ax7 = plt.subplots(constrained_layout=True)
        ax7.plot(ts, delta_prod, marker = 'o', color='royalblue', mfc='white', label=rf'N$_{{steps}}$ = {Nsteps[j]}')
        # ax7.plot(ts, inc_prod, marker = 'o', color='deeppink', mfc='white', label=rf'N_{{steps}} = {Nsteps[j]}')
        ax7.axhline(y=hbar/2, color='black', linestyle = '--', label=r'$\frac{\hbar}{2}$')
        ax7.set_xlabel(r't [s]', fontsize=fs)
        ax7.set_ylabel(r"$\langle \Delta x \rangle(t) \langle \Delta p \rangle(t)$", fontsize=fs)
        ax7.tick_params(axis="both", labelsize=fs)
        ax7.legend(fontsize=fs, loc='lower left')
        ax7.grid(True)


plt.show()




