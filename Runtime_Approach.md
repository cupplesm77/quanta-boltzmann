The main slowdown is not the array itself. It is the 100 million-iteration Python loop, plus two 
relatively expensive np.random.choice() calls in many iterations.

For this sequential simulation, the biggest improvement will come from compiling the loop with Numba. 
Numba’s @njit is specifically intended to compile numerical loops into machine code without running 
each iteration through the Python interpreter.

Recommended optimized version
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from numba import njit
from time import perf_counter
the_seed = 78781


@njit
def run_simulation(boltz_array, num_iterations, seed=the_seed):
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


length_boltz_array = 100_000
num_iterations = 100 * 1_000_000

# These values are counts, so use an integer dtype from the beginning
boltz_array = np.ones(length_boltz_array, dtype=np.int64)

initial_total = boltz_array.sum()

start_time = perf_counter()

boltz_array = run_simulation(
    boltz_array,
    num_iterations,
    seed=the_seed
)

elapsed_time = perf_counter() - start_time

print(f"Runtime: {elapsed_time:.3f} seconds")
print(f"Initial total: {initial_total}")
print(f"Final total:   {boltz_array.sum()}")
print(f"Minimum value: {boltz_array.min()}")
print(f"Maximum value: {boltz_array.max()}")

Count the values efficiently

# Since all the values are nonnegative integers, np.bincount() is a natural replacement for 
np.unique(..., return_counts=True):

counts = np.bincount(boltz_array)

# Integer population values represented by the positions in counts
values = np.arange(counts.size)

# Optional: remove trailing or internal zero-count values from the plot
occupied = counts > 0
values = values[occupied]
counts = counts[occupied]


Then plot:

fig, ax = plt.subplots()

ax.plot(values, counts, marker="o")

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_xlabel("Number of particles at a site")
ax.set_ylabel("Number of sites")
ax.set_title("Particle-count distribution")

plt.show()

Why this should be much faster
1. Replace np.random.choice(idx)

You currently have:

idx = range(len(botz_array))
random_index1 = np.random.choice(idx)


But you only need a uniformly distributed integer index. The direct operation is:

random_index1 = np.random.randint(length_boltz_array)


Outside Numba, the modern NumPy equivalent is:

rng = np.random.default_rng()
random_index1 = rng.integers(length_boltz_array)


Generator.integers(high) draws an integer from zero, inclusive, to high, exclusive.

2. Eliminate operations that do nothing

These statements can simply be removed:

botz_array[random_index2] += 0
botz_array[random_index1] += 0


They do not alter the array, but they still require indexing and assignment.

Your logic can be reduced to:

if boltz_array[index1] > 0:
    index2 = ...
    if index2 != index1:
        boltz_array[index1] -= 1
        boltz_array[index2] += 1

3. Use integer storage from the beginning

You currently have:

botz_array = np.ones(length_boltz_array)


By default, that normally creates floating-point values. Since these are particle counts, create integers 
directly:

boltz_array = np.ones(
    length_boltz_array,
    dtype=np.int64
)


Then you no longer need:

unq = unique.astype(int)

4. Compile the loop

Normal NumPy speedups often come from replacing loops with vectorized array operations. This particular 
loop is difficult to vectorize correctly because iteration j+1j+1 depends on changes made during iteration jj.

For example, generating 100 million index pairs and applying them all simultaneously would not automatically 
reproduce the same process. Repeated indices, empty source sites, and the ordering of transfers all matter.

This is exactly the kind of case where retaining the loop but compiling it with Numba is appropriate. Numba’s 
documentation notes that loops work well under @njit, while the same loops can be much slower in ordinary Python.

Important first-run detail

The first call:

run_simulation(...)


includes Numba’s compilation overhead. Subsequent calls using the same argument types reuse the compiled version.

For proper timing, warm it up first with a tiny simulation:

warmup = np.ones(10, dtype=np.int64)
run_simulation(warmup, 10)


Then time the real run:

start_time = perf_counter()
run_simulation(boltz_array, num_iterations)
elapsed_time = perf_counter() - start_time

A useful correctness check

Every successful move subtracts one from one site and adds one to another. Therefore, the total should 
remain constant:

assert boltz_array.sum() == length_boltz_array
assert boltz_array.min() >= 0


Those two assertions verify:

Particle number is conserved.
No site has gone below zero.

One small naming correction would also help readability: change botz_array to boltz_array. The optimized 
code preserves the sequential dynamics while removing unnecessary work and compiling the expensive loop.