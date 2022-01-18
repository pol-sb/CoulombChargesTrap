import gc
import os
import time

import matplotlib.pyplot as plt
import numpy as np

# TODO: Add docstrings to all of the functions!

# TODO: The monte carlo function should be optimized so that not all
# interactions are computed at each iteration, only the interactions of the
# excited particle with the rest. This is correct, right? Check.

# TODO: Simplify the result treatment function.


def print_color(text: str, color: str):
    d_col = {
        "blue": "\033[1;34m",
        "green": "\033[1;32m",
        "red": "\033[1;31m",
        "yellow": "\033[1;33m",
        "normal": "\033[0m",
    }
    print(f"{d_col[color]}{text}{d_col['normal']}")


def write_def_xyz(
    filename: str,
    pos_arr: np.ndarray,
    n_part: int,
    traj: bool = True,
):
    with open(f"{filename}", "w+") as f:
        if traj:
            for pos_v in pos_arr:
                f.write(str(n_part) + "\n")
                f.write("\n")
                for pos in pos_v:
                    f.write(f"A\t{pos[0]}\t{pos[1]}\t{pos[2]}\n")
        else:
            f.write(str(n_part) + "\n")
            f.write("\n")
            for pos in pos_arr:
                f.write(f"A\t{pos[0]}\t{pos[1]}\t{pos[2]}\n")


# n>1 initial conditions should be generated to be reasonably certain that
# all possible local minima have been found.
# The literature uses 50.
def initialize_system(
    n_part: int,
    n_confs: int,
    L: float,
    rng: np.random.Generator,
):
    conf_list = []
    for conf in range(n_confs):
        init_coords = L * (2 * rng.random([n_part, 3]) - 1)
        conf_list.append(init_coords)
    return conf_list

# TODO: Optimize this function!
def total_pot_ener(cor: np.ndarray, **params):

    alpha = params["alpha"]
    d = params["d"]
    q = params["q"]
    L = params["L"]

    sum_term = 0
    for i in range(cor.shape[0]):
        sum_term_j = 0

        xi, yi, zi = cor[i, 0], cor[i, 1], cor[i, 2]
        # term_d = 1.0 / (2.0 * (d ** 3))
        # term_coord = (zi ** 2) + alpha * ((xi ** 2) + (yi ** 2))
        term = (xi ** 2 + yi ** 2 + zi ** 2) / 2

        # sum_term += term_d * term_coord
        sum_term += term
        for j in range(i):
            dist = np.linalg.norm((cor[i, :] - cor[j, :]))
            # dist = 1 / distancesq(L, cor[i, :], cor[j, :])
            # modu = 1 / np.sum(np.abs((cor[i, :] - cor[j, :])))
            sum_term_j += 1 / dist

        sum_term += sum_term_j
        # print('sum_term: ', sum_term)

    tot_ene = (q ** 2) * sum_term

    return tot_ene


# The Monte Carlo algorithm, used principally for N < 10 and a != 1,
# moves ions randomly within a step size which depends on the T and keeps the
# move only if the energy of the configuration is reduced or a random check in a
# Maxwell-Bolzmann distribution is passed.
def monte_carlo(rng: np.random.Generator, **params):

    n_MCS = int(params["n_MCS"])
    coef_T = params["coef_T"]

    res_dict = {}

    for N in params["N_list"]:

        coord_arr = initialize_system(
            N,
            params["confs"],
            params["L"],
            rng,
        )

        final_ener_list = []
        final_T_list = []
        final_pos_list = []
        radius_list = []

        print_color(
            f"\n[{time.strftime('%H:%M:%S')}] - Working with N = {N} particles",
            "green",
        )

        for c_ind, conf in enumerate(coord_arr):
            T = params["T"]
            time_a = time.time()
            Nacc = 0

            print_color(
                f"\n[{time.strftime('%H:%M:%S')}] - Configuration"
                f" {c_ind+1}/{len(coord_arr)}:",
                "blue",
            )

            print(f"Initial energy:\n{total_pot_ener(conf, **params)}\n")
            print("Progress:")

            conf_traj = []
            energ_list = []
            T_list = []

            syst_shape = conf.shape

            T_list.append(T)
            initial_E = total_pot_ener(conf, **params)
            energ_list.append(initial_E)

            for mcs in range(n_MCS):

                # Printing progress and energy every n iterations.
                if mcs % 1000 == 0:
                    print(f"{mcs}/{n_MCS} - Energy: {initial_E:.7f}", end="\r")

                # Preparing the random perturbation for each particle.
                # A single perturbation is applied to a single random particle,
                # so the MC acceptance is higher.
                perturb = np.zeros(syst_shape)
                ind_col = rng.integers(0, N)
                perturb[ind_col, :] = T * (2 * rng.random(3) - 1)

                # Applying the perturbation
                conf_next = conf + perturb

                # Computing the energy after the perturbation
                perturb_E = total_pot_ener(conf_next, **params)

                # Difference in energy between initial E and perturbed E.
                delta_E = perturb_E - initial_E

                # Accept the new geometry if the energy is reduced
                if delta_E < 0:
                    Nacc += 1
                    initial_E = perturb_E
                    conf = conf_next

                # If the energy is not reduced, use the probability distrib.
                # to check if the change is accepted anyway.
                else:
                    rand = rng.random()
                    if np.exp(-delta_E / T) > rand:
                        Nacc += 1
                        initial_E = perturb_E
                        conf = conf_next

                # The coefficient coef_T reduces the T each iteration in a way
                # that will result in reaching the final T in the last iteration
                # of the execution.
                T *= coef_T

                # Writing each configuration, T and E to a list, which will
                # be plotted or exported later to a '.xyz' file for
                # representation.
                energ_list.append(perturb_E)
                T_list.append(T)
                conf_traj.append(conf)

            Pacc = (Nacc / n_MCS) * 100

            time_b = time.time()
            tot_t = time_b - time_a
            print_color("\n\nResults:", "blue")

            radius = 0
            for part in conf:
                radius += part[0] ** 2 + part[1] ** 2 + part[2] ** 2
            radius /= syst_shape[0]
            radius = np.sqrt(radius)

            print("Radius:", radius)
            print("Final energy:", initial_E)
            print(f"N_accept: {Nacc}\t\t P_accept: {Pacc:.3f}%")
            print("Final T:", T)
            print_color(
                f"[{time.strftime('%H:%M:%S')}] - Done. Time elapsed:"
                f" {tot_t:.2f}s.",
                "blue",
            )

            radius_list.append(radius)
            final_ener_list.append(energ_list)
            final_pos_list.append(conf_traj)
            final_T_list.append(T_list)
            gc.collect()

        res_dict[f"{N}"] = {
            "Energies": final_ener_list,
            "Temperatures": final_T_list,
            "Trajectory": final_pos_list,
            "Radii": radius_list,
        }

        # Freeing memory?
        gc.collect()

        # Clearing lists.
        # TODO: Does this serve any purpose? I think it does not clear
        # any memory.
        conf_traj = []
        energ_list = []
        T_list = []
        final_ener_list = []
        final_T_list = []
        final_pos_list = []
        radius_list = []

    return res_dict


def result_treatment(res_dict, sys_args):

    plt.style.use("seaborn-poster")
    plt.switch_backend("agg")

    print_color("\n\nCalculations Done.", "green")
    print_color("Saving results...\n", "green")

    c_time = time.strftime("%d-%H_%M_%S")
    base_path = f"iontrap_{c_time}"
    os.mkdir(base_path)

    for N in res_dict.keys():

        path_N = f"{base_path}/N_{N}"
        os.mkdir(path_N)
        for c_ind in range(len(res_dict[f"{N}"]["Energies"])):

            final_conf = res_dict[f"{N}"]["Trajectory"][c_ind]
            E_list = res_dict[f"{N}"]["Energies"][c_ind]
            T_list = res_dict[f"{N}"]["Temperatures"][c_ind]
            radius = res_dict[f"{N}"]["Radii"][c_ind]
            n_MCS = int(sys_args["n_MCS"])

            # Saving the final configuration (or the entire system trajectory)
            # to a .xyz file.
            write_def_xyz(
                f"{path_N}/conf_{c_ind+1}_traj.xyz",
                final_conf,
                N,
                traj=True,
            )

            with open(f"{path_N}/conf_{c_ind+1}_results.resmc", "w+") as f:
                f.write(
                    f"### RESULTS FOR N = {N}, CONFIGURATION {c_ind+1} ###\n\n"
                )
                f.write(f"Total iterations: {n_MCS}\n")
                f.write(f"Final Energy: {E_list[-1]}\n")
                f.write(f"Radius: {radius}\n")
                f.write(f"Final T: {T_list[-1]}\n")

            # Plotting the energies over time.
            fig = plt.figure(num=1, clear=True)
            ax = fig.add_subplot()
            ax.plot(range(len(E_list)), E_list)
            ax.set_xlabel("Iterations")
            ax.set_ylabel("Energy")
            ax.set_title(f"Evolution of the energy for N = {N}")
            fig.savefig(f"{path_N}/conf_{c_ind+1}_E_v_iter.svg", dpi=800)
            fig.savefig(
                f"{path_N}/conf_{c_ind+1}_E_v_iter_nobg.svg",
                transparent=True,
            )

            fig = plt.figure(num=1, clear=True)
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.plot(range(len(E_list)), E_list)
            ax1.set_xlabel("Iterations")
            ax1.set_ylabel("Energy")
            ax1.set_title(f"Evolution of the energy for N = {N}")

            ax2 = fig.add_subplot(2, 1, 2)
            ax2.plot(
                range(int(n_MCS * 0.60), n_MCS + 1),
                E_list[int((n_MCS * 0.60)) :],
            )
            ax2.set_title(
                f"Evolution of the energy for N = {N}, after"
                f" {int(n_MCS * 0.60)} iter."
            )
            ax2.set_xlabel("Iterations")
            ax2.set_ylabel("Energy")
            fig.tight_layout()
            fig.savefig(f"{path_N}/conf_{c_ind+1}_E_v_iter_doble.svg", dpi=800)
            fig.savefig(
                f"{path_N}/conf_{c_ind+1}_E_v_iter_doble_nobg.svg",
                transparent=True,
            )

            fig = plt.figure(num=1, clear=True)
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.plot(T_list, E_list)
            ax1.set_xlim(max(T_list), 0)
            ax1.set_xlabel("T")
            ax1.set_ylabel("Energy")
            ax1.set_title(f"Temperature effect on the Energy (N = {N})")

            ax2 = fig.add_subplot(2, 1, 2)
            ax2.plot(
                T_list[int((n_MCS * 0.60)) :],
                E_list[int((n_MCS * 0.60)) :],
            )
            ax2.set_xlim(max(T_list[int((n_MCS * 0.60)) :]), 0)
            ax2.set_title(
                f"Temperature effect on the Energy (N = {N}), after"
                f" {int(n_MCS * 0.60)} iter."
            )
            ax2.set_xlabel("T")
            ax2.set_ylabel("Energy")
            fig.tight_layout()
            fig.savefig(
                f"{path_N}/conf_{c_ind+1}_E_v_T.svg",
                transparent=False,
            )
            fig.savefig(
                f"{path_N}/conf_{c_ind+1}_E_v_T_nobg.svg",
                transparent=True,
            )

            # Garbage Collection
            gc.collect()

        en_list = res_dict[f"{N}"]["Energies"]
        t_list = res_dict[f"{N}"]["Temperatures"]

        last_E_list = np.array([c_list[-1] for c_list in en_list])
        min_E = np.min(last_E_list)
        best_c = np.where(last_E_list == min_E)[0][0]

        with open(f"{path_N}/best_run_results.resmc", "w+") as f:
            f.write(f"### BEST CONFIGURATION FOR N = {N} ###\n\n")
            f.write(f"CONFIGURATION {best_c+1}\n")
            f.write(f"Total iterations: {n_MCS}\n")
            f.write(f"Final Energy: {en_list[best_c][-1]}\n")
            f.write(f"Radius: {radius}\n")
            f.write(f"Final T: {t_list[best_c][-1]}\n")

    print_color(f"Results stored in '{base_path}'.\n", "green")


class SimulateIonsTrap:
    def __init__(self, sys_arg):
        rng = np.random.default_rng()

        res_dict = monte_carlo(rng, **sys_arg)
        result_treatment(res_dict, sys_arg)
        gc.collect()


if __name__ == "__main__":
    _red = "\033[1;31m"
    _nc = "\033[0m"
    print(
        f"\n{_red}[!!] ERROR:{_nc} This module should only be imported and not"
        " run directly.\n"
    )
