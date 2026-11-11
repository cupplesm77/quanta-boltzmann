# statistical nature of the boltzmann distribution

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from numba import njit
from time import perf_counter
boltz_seed = 78781

@njit
def run_simulation(boltz_array, num_iterations, seed=boltz_seed):
    np.random.seed(seed)

    array_length = boltz_array.size

    for _ in range(num_iterations):

        # Faster and clearer than np.random.choice(range(...))
        index1 = np.random.randint(array_length)

        # A transfer is possible only if this site is occupied
        if boltz_array[index1] > 0:

            index2 = np.random.randint(array_length)

            # Do nothing if both indices refer to the same site
            if index2 != index1:
                boltz_array[index1] -= 1
                boltz_array[index2] += 1

    return boltz_array

# boltzmann array initialization routine
def initialize_boltzmann(Q, M):
    q, r = divmod(Q, M)

    # Give every lattice point q quanta
    boltzmann = np.full(M, q, dtype=np.int64)

    # Spread the remaining r quanta approximately uniformly
    if r > 0:
        indices = (np.arange(r) * M) // r
        boltzmann[indices] += 1

    return boltzmann


M = 100_000  # number of lattice sites
Q = 100_000  # total number of energy quanta

million = 1_000_000
num_iterations = 10 * million

# values are counts of energy quanta at each lattice site, so use an integer dtype
# boltz_array = np.ones(M, dtype=np.int64)

# Construct an initialization routine that randomly disperses the Q quanta over the M lattice points
boltz_array0 = initialize_boltzmann(Q, M)
assert len(boltz_array0) == M

initial_total = boltz_array0.sum()
initial_size = len(boltz_array0)

# For proper timing, warm simulation up first with a tiny run:
warmup = np.ones(10, dtype=np.int64)
_ = run_simulation(warmup, 10)


# Primary Simulation Loop: Run Multiple Simulations
sim_runs = 1
run_times = []
# inc = 100_000
inc = 0

for run in range(sim_runs):

    # time the simulation run
    start_time = perf_counter()

    boltz_array = boltz_array0.copy()
    run_seed = boltz_seed  + run

    print("")
    print("CHECK: Beginning run:", run)
    print("CHECK: Initial first 10 sites:", boltz_array[:10])
    print("CHECK: Initial total Q:", boltz_array.sum())

    # add increments of interation to each of the simulation runs

    boltz_array = run_simulation(
        boltz_array,
        num_iterations,
        seed=run_seed
    )
    num_iterations += inc

    print("")
    print("CHECK: Final first 10 sites:", boltz_array[:10])
    print("CHECK: Final total Q:", boltz_array.sum())
    print("")

    # Every successful move subtracts one from one site and adds one to another.
    # Therefore, the total quanta should remain constant:

    # no lattice points were lost
    assert len(boltz_array) == M
    # conservation of quanta (weak)
    assert boltz_array.sum() == Q
    # no negative lattice values
    assert boltz_array.min() >= 0

    elapsed_time = perf_counter() - start_time
    run_times.append(elapsed_time)

    print(f"Runtime: {elapsed_time:.3f} seconds")
    print(f"Initial total: {initial_total}")
    print(f"Final total:   {boltz_array.sum()}")
    print(f"Minimum value: {boltz_array.min()}")
    print(f"Maximum value: {boltz_array.max()}")

    # Since all the values are nonnegative integers, np.bincount() works nicely
    counts = np.bincount(boltz_array)
    print("counts = ", counts)

    # Integer population values represented by the positions in counts
    values = np.arange(counts.size)
    print("values = ", values)
    print("")

    # strong conservation of quanta
    assert (counts * values).sum() == Q

    # Remove trailing or internal zero-count values from the plot
    occupied = counts > 0
    values = values[occupied]
    counts = counts[occupied]

    fig, ax = plt.subplots()
    ax.plot(values, counts, marker="o", linestyle="None")
    # use a logarithmic scale for the y-axis
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Number (n) of particles at a site")
    ax.set_ylabel("Number (N) of sites with n particles")
    ax.set_title("Particle-count distribution")

    plt.show()


