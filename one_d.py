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
Q = 200_000  # total number of energy quanta

num_iterations = 300 * 1_000_000

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
sim_runs = 3
run_times = []
inc = 10
for run in range(sim_runs):

    # time the simulation run
    start_time = perf_counter()

    # add increments of interation to each of the simulation runs
    num_iterations += inc
    boltz_array = run_simulation(
        boltz_array0,
        num_iterations,
        seed=boltz_seed
    )

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


# Slow Approach to Problem:
if False:
    M = 100000
    num_iterations = 100 * 1e6


    botz_array = np.ones(M)
    print("size of boltz array: ", botz_array.size)
    print(botz_array)

    # build a routine that will generate a boltzmann distribution
    print(botz_array.ndim)
    print(botz_array[1:4])

    # for x in botz_array:
    #     print(x)

    idx = range(len(botz_array))
    print("idx: ", idx)
    # print("")

    # iterate the boltzmann array
    test_arr = []
    n=0
    for j in range(num_iterations):

        random_index1 =  np.random.choice(idx)
        # print(f"random_index1: {random_index1}")
        x = botz_array[random_index1]
        # print(f"random element of array: {x}")
        if x != 0:
            #print("x is not zero")
            # choose another random index and add x to the value
            random_index2 =  np.random.choice(idx)

            if random_index2 != random_index1:
                botz_array[random_index2] += 1
                botz_array[random_index1] -= 1
                # test_arr.append(random_index1)
                # n = n + 1
            else:
                # don't do anything if we choose the indexes of the lattice point equal to the idx value
                botz_array[random_index2] += 0
                # continue
        else:
            # if there is nothing in the lattice point, do nothing
            botz_array[random_index1] += 0
            # continue


    print(botz_array)
    unique, counts = np.unique(botz_array, return_counts=True)
    unq = unique.astype(int)
    print("unique = ", unq)
    print("counts = ", counts)
    print("test_arr size = ", len(test_arr))
    print("n = ", n)

    # plot the boltzmann distribution
    # add proper axis labels
    plt.plot(unq, counts)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel("unique values")
    plt.ylabel("counts")
    plt.show()


