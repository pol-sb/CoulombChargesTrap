# Following: Rafac, R., Schiffer, J. P., Hangst, J. S., Dubin, D. H., &
# Wales, D. J. (1991). Stable configurations of confined cold ionic systems.
# In Proceedings of the National Academy of Sciences
# (Vol. 88, Issue 2, pp. 483–486).

import mc_utils as mcu

# Dictionary with parameters for the simulation
# N_list: List of number of particles.
# n_MCS: Number of iterations for the Monte Carlo method
# L: Lattice size -> L
# q: Ionic charge
# alpha: Symmetry of the external harmonic potentia
# m: Mass
# omega: Natural frequency
ar_dict = {
    "n_MCS": 1e5,
    "L": 2,
    "q": 1,
    "m": 1,
    "alpha": 0.5,
    "omega": 1,
    "T": 10,
    "T_f": 10e-6,
    "confs": 2,
    "plot": False,
}

N_list = [4]

# Unit length (d)
d = ((ar_dict["q"] ** 2) / ar_dict["m"] * (ar_dict["omega"] ** 2)) ** (1 / 3)
ar_dict["d"] = d

# Computing the coefficient for the T
coef_T = (ar_dict["T_f"] / ar_dict["T"]) ** (1 / ar_dict["n_MCS"])
ar_dict["coef_T"] = coef_T

# Runs the simulation with the given parameters
for N in N_list:
    ar_dict["N_list"] = [N]
    mcu.SimulateIonsTrap(ar_dict)
