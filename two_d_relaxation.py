# Statistical nature of the Boltzmann distribution on a 2-D lattice
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from numba import njit
from time import perf_counter

boltz_seed = 78781


@njit
def run_simulation_with_snapshots(
        boltz_array,
        snapshot_iterations,
        seed=boltz_seed,
):
    """
    Simulate random energy quanta exchange on a 2D lattice with snapshot logging.

    This function models the relaxation of an isolated system of discrete
    harmonic oscillators (an Einstein solid) toward statistical equilibrium.
    At each iteration, energy is conserved through pairwise random exchange:
    a donor site is chosen uniformly at random and, if it possesses at least
    one quantum (E >= 1), transfers one quantum to a randomly chosen receiver
    site. Over many iterations, random exchanges maximize system entropy and
    drive the single-site occupancy distribution toward the Boltzmann
    (geometric) distribution.

    Parameters
    ----------
    boltz_array : numpy.ndarray
        2D integer array of shape (Lx, Ly) representing the initial lattice
        microstate, where each element stores the number of energy quanta at
        that site. This array is modified in-place during the simulation.
    snapshot_iterations : numpy.ndarray
        1D integer array of strictly increasing iteration counts (1-based)
        at which to record snapshots of the 2D lattice configuration.
    seed : int, optional
        Seed for the random number generator to ensure reproducible
        simulation trajectories. Default is `boltz_seed`.

    Returns
    -------
    boltz_array : numpy.ndarray
        The modified 2D lattice array representing the final microstate of the
        system after the last requested iteration.
    snapshots : numpy.ndarray
        3D integer array of shape (N, Lx, Ly), where N is the number of
        requested checkpoints, containing copies of the 2D lattice state
        captured at each iteration in `snapshot_iterations`.

    Notes
    -----
    - Physical Conservation: Total quanta Q = sum(boltz_array) and lattice
      size M = Lx * Ly are strictly conserved throughout the simulation.
    - Global Exchange: Energy transfers occur between any two sites across
      the entire lattice (mean-field interaction), rather than strictly nearest
      neighbors, enabling rapid mixing and relaxation toward equilibrium.
    """
    np.random.seed(seed)

    # flat is a view of the same data, not a separate physical system.
    # This preserves the original global-transfer dynamics.
    flat = boltz_array.ravel()
    array_length = flat.size

    snapshots = np.empty(
        (snapshot_iterations.size, boltz_array.shape[0], boltz_array.shape[1]),
        dtype=boltz_array.dtype,
    )

    snapshot_index = 0
    final_iteration = snapshot_iterations[-1]

    for iteration in range(1, final_iteration + 1):
        index1 = np.random.randint(array_length)

        # A transfer is possible only if the source site is occupied.
        if flat[index1] > 0:
            index2 = np.random.randint(array_length)

            # Do nothing if both indices refer to the same site.
            if index2 != index1:
                flat[index1] -= 1
                flat[index2] += 1

        if iteration == snapshot_iterations[snapshot_index]:
            snapshots[snapshot_index, :, :] = boltz_array
            snapshot_index += 1
            if snapshot_index == snapshot_iterations.size:
                break

    return boltz_array, snapshots


def initialize_boltzmann(Q, Lx, Ly):
    """
    Distribute energy quanta across a 2D lattice in a uniform, low-entropy state.

    Calculates the base number of quanta per site q = Q // M and remainder
    r = Q % M for a lattice of total size M = Lx * Ly. Every site receives at
    least q quanta, and the remaining r quanta are dispersed as evenly as
    possible across the lattice sites. This creates a deterministic, minimally
    disordered initial configuration for studying relaxation toward statistical
    equilibrium.

    Parameters
    ----------
    Q : int
        Total number of discrete energy quanta to distribute across the lattice.
    Lx : int
        Number of lattice rows (dimension along the first axis).
    Ly : int
        Number of lattice columns (dimension along the second axis).

    Returns
    -------
    boltzmann : numpy.ndarray
        2D integer array of shape (Lx, Ly) and dtype `numpy.int64` representing
        the initial lattice microstate with sum(boltzmann) == Q.

    Notes
    -----
    - Statistical Mechanics Context: Initializing the system with a nearly
      uniform distribution (minimal variance in local energy) establishes a
      low-entropy microstate. As the simulation runs, random exchanges increase
      entropy until the single-site occupancies reach the Boltzmann (geometric)
      distribution.
    - Uniform Remainder Spacing: When Q is not an exact multiple of M, the r
      excess quanta are assigned to flattened 1D indices `(k * M) // r` for
      `k = 0, 1, ..., r - 1`, ensuring balanced spatial dispersal across the
      lattice.
    """
    # Count the number of oscillators, M = Lx * Ly .
    M = Lx * Ly

    # Divide the Q energy quanta among the M oscillators as evenly as possible.
    q, r = divmod(Q, M)

    # Give every oscillator q quanta using np.full().
    boltzmann = np.full((Lx, Ly), q, dtype=np.int64)

    # Spread any remaining quanta approximately uniformly using a flat view.
    if r > 0:
        # If there are r quanta left over, use ravel() to treat the 2-D lattice as a simple 1-D sequence of sites.
        flat = boltzmann.ravel()
        # Choose r roughly evenly spaced sites and give each one additional quantum.
        indices = (np.arange(r) * M) // r
        flat[indices] += 1

    assert boltzmann.sum() == Q

    return boltzmann


def equillibrium_dist(boltz_array, rho):
    """Return total-variation distance from geometric equilibrium."""
    counts = np.bincount(boltz_array.ravel())
    p_measured = counts / boltz_array.size
    n = np.arange(counts.size)
    ratio = rho / (1.0 + rho)
    p_equilibrium = (1.0 / (1.0 + rho)) * ratio**n

    # p_measured is zero above the largest observed occupation.
    # Add the exact geometric probability remaining in that tail.
    equilibrium_tail = ratio**counts.size
    return 0.5 * (np.abs(p_measured - p_equilibrium).sum() + equilibrium_tail)


@njit
def run_relaxation_history(boltz_array, total_sweeps, seed, rho):
    """Run one independent experiment and record D(t) after every sweep."""
    np.random.seed(seed)
    flat = boltz_array.ravel()
    M = flat.size
    distances = np.empty(total_sweeps + 1, dtype=np.float64)

    # D(0): exact initial condition.
    # Inline the distance calculation so this time-critical loop remains JIT compiled.
    max_n = flat.max()
    counts = np.zeros(max_n + 1, dtype=np.int64)
    for value in flat:
        counts[value] += 1
    ratio = rho / (1.0 + rho)
    prefactor = 1.0 / (1.0 + rho)
    distance = 0.0
    for n in range(counts.size):
        p_measured = counts[n] / M
        p_equilibrium = prefactor * ratio**n
        distance += abs(p_measured - p_equilibrium)
    distance += ratio**counts.size
    distances[0] = 0.5 * distance

    for sweep in range(1, total_sweeps + 1):
        # One sweep = M attempted updates.
        for _ in range(M):
            index1 = np.random.randint(M)
            if flat[index1] > 0:
                index2 = np.random.randint(M)
                if index2 != index1:
                    flat[index1] -= 1
                    flat[index2] += 1

        max_n = flat.max()
        counts = np.zeros(max_n + 1, dtype=np.int64)
        for value in flat:
            counts[value] += 1
        distance = 0.0
        for n in range(counts.size):
            p_measured = counts[n] / M
            p_equilibrium = prefactor * ratio**n
            distance += abs(p_measured - p_equilibrium)
        distance += ratio**counts.size
        distances[sweep] = 0.5 * distance

    return boltz_array, distances


# Lattice geometry and conserved total
Lx = 100
Ly = 100
M = Lx * Ly
Q = M                         # density rho = Q/M = 1
rho = Q / M

# Each value is now an independent relaxation experiment.
# num_sweeps = [1, 10, 100, 500, 1000, 2000, 5000]
num_sweeps = [1, 10, 20, 30, 50, 100, 1000]

boltz_array0 = initialize_boltzmann(Q, Lx, Ly)
assert boltz_array0.shape == (Lx, Ly)
assert boltz_array0.size == M
initial_total = boltz_array0.sum()

# Warm up the JIT-compiled relaxation routine before timing experiments.
warmup = np.ones((2, 5), dtype=np.int64)
_ = run_relaxation_history(warmup, 1, boltz_seed, 1.0)

final_snapshots = []
distance_histories = []
run_times = []

for case_index, total_sweeps in enumerate(num_sweeps):
    # Every snapshot case starts independently from the same pristine state.
    boltz_array = boltz_array0.copy()
    run_seed = boltz_seed + case_index

    print()
    print(f"CHECK: Beginning independent {total_sweeps:,}-sweep case")
    start_time = perf_counter()

    boltz_array, distance_data = run_relaxation_history(
        boltz_array,
        total_sweeps,
        run_seed,
        rho,
    )

    elapsed_time = perf_counter() - start_time
    run_times.append(elapsed_time)
    final_snapshots.append(boltz_array.copy())
    distance_histories.append(distance_data)

    assert boltz_array.shape == (Lx, Ly)
    assert boltz_array.size == M
    assert boltz_array.sum() == Q
    assert boltz_array.min() >= 0
    assert distance_data.size == total_sweeps + 1
    assert np.isclose(distance_data[0], equillibrium_dist(boltz_array0, rho))

    print(f"Runtime: {elapsed_time:.3f} seconds")
    print(f"D(0):                 {distance_data[0]:.6f}")
    print(f"D({total_sweeps:,} sweeps): {distance_data[-1]:.6f}")


# Graphics
# Include the pristine initial state before the independently evolved cases.
snapshots = np.concatenate(
    (
        boltz_array0[np.newaxis, :, :],
        np.stack(final_snapshots),
    ),
    axis=0,
)

snapshot_labels = ["Initial"] + [f"{n:,} sweeps" for n in num_sweeps]

# Plot the initial state followed by the final spatial state
# from every independent experiment.
vmax = snapshots.max()
fig, axes = plt.subplots(
    1,
    snapshots.shape[0],
    figsize=(20, 4),
    constrained_layout=True,
)
for ax, state, label in zip(axes, snapshots, snapshot_labels):
    image = ax.imshow(
        state,
        origin="lower",
        cmap="viridis",
        vmin=0,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_title(label)
    ax.set_xlabel("y site")
    ax.set_ylabel("x site")
colorbar = fig.colorbar(image, ax=axes, shrink=0.82)
colorbar.set_label("Quanta at site n(x, y)")
fig.suptitle(f"Independent 2-D relaxation experiments: {Lx} x {Ly}, Q = {Q:,}")

# One independent D(t) plot for every snapshot case.
fig2, axes2 = plt.subplots(2, 4, figsize=(15, 8), constrained_layout=True)
axes2 = axes2.ravel()
for ax, total_sweeps, distance_data in zip(axes2, num_sweeps, distance_histories):
    sweep_times = np.arange(total_sweeps + 1)
    ax.plot(sweep_times, distance_data, marker="o" if total_sweeps <= 10 else None)
    ax.set_xlabel("Time (sweeps)")
    ax.set_ylabel("D(t)")
    ax.set_title(f"run: {total_sweeps:,} sweeps")
    ax.grid(True, alpha=0.3)

# Turn off only unused subplot slots
for ax in axes2[len(num_sweeps):]:
    ax.axis("off")

fig2.suptitle("Evolution of relaxation toward equilibrium")

# Retain the final particle-count distribution, using the longest independent run.
boltz_array = snapshots[-1]
counts = np.bincount(boltz_array.ravel())
values = np.arange(counts.size)
sigma_poisson = np.sqrt(counts)
assert (counts * values).sum() == Q
assert counts.sum() == M
occupied = counts > 0

fig3, ax3 = plt.subplots()
ax3.errorbar(
    values[occupied],
    counts[occupied],
    yerr=sigma_poisson[occupied],
    marker="o",
    linestyle="None",
    capsize=3,
)
ax3.set_yscale("log")
ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
ax3.set_xlabel("Number (n) of particles at a site")
ax3.set_ylabel("Number (N) of sites with n particles")
ax3.set_title(f"Final particle-count distribution ({num_sweeps[-1]:,}-sweep run)")

plt.show()
