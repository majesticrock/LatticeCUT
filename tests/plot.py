import numpy as np
import matplotlib.pyplot as plt
import mrock_centralized_scripts.path_appender as __ap
__ap.append()

import mrock_centralized_scripts.FullDiagPurger as fdp
from get_data import load_panda, lattice_cut_params
import sys

def load_data(name):
    """ Returns the data for the LatticeCUT test with the name 'name'.
    standard: SC lattice g=2.5 (for secondary excitations)
    enhanced: BCC with g=1.5 U=0.2 (for enhanced mode)
    """
    
    if name == "standard":
        return load_panda("lattice_cut", "test/sc", "full_diagonalization.json.gz",
                    **lattice_cut_params(N=2000, g=2.5, U=0., E_F=0, omega_D=0.02))
    elif name == "enhanced":
        return load_panda("lattice_cut", "test/bcc", "full_diagonalization.json.gz",
                    **lattice_cut_params(N=2000, g=1.5, U=0.2, E_F=-0.5, omega_D=0.02))
    else:
        raise ValueError("Continuum test: Invalid index")

def create_plot(name):
    main_df = load_data(name)
    
    if name=="standard":
        fig, axes = plt.subplots(nrows=3, sharex=True)
        fig.subplots_adjust(hspace=0)
        axes[0].set_ylabel("Higgs")
        axes[1].set_ylabel("Occupation")
        axes[2].set_ylabel("Phase")
        axes[-1].set_xlabel(r"$\varepsilon - \mu$")
        axes[0].set_xlim(-0.25, 0.25)
        purger = fdp.FullDiagPurger(main_df, np.linspace(-1, 1, 2000) - main_df["chemical_potential"])
    
        purger.plot_amplitude(axes[:2], combined_norm=True)
        purger.plot_phase(axes[2], label="Result")
    elif name=="enhanced":
        fig, ax = plt.subplots()
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(r"$A (\omega)$")
        
        omegas = np.linspace(0., 0.4, 5000, dtype=complex) + 5e-3j
        A_phase = -np.imag(np.array([np.sum(main_df["phase.weights"] / (z**2 - main_df["phase.eigenvalues"]**2)) 
                                for z in omegas])) / np.pi
        A_higgs = -np.imag(np.array([np.sum(main_df["amplitude.weights"] / (z**2 - main_df["amplitude.eigenvalues"]**2)) 
                                for z in omegas])) / np.pi
    
        ax.plot(omegas.real, A_higgs, label="Higgs")
        ax.plot(omegas.real, A_phase, label="Phase")
        ax.set_ylim(0, 0.7)
        ax.legend()
    
    fig.suptitle(f"LatticeCUT: {name}")
    
    plt.show()

if len(sys.argv) > 1:
    create_plot(sys.argv[1])
else:
    print("Please provide the name of the test you would like to plot.")