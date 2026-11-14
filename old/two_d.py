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

            # FUNCTION HERE  (or Choose a better placement).

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


# Lattice geometry and conserved total
Lx = 100
Ly = 100
M = Lx * Ly
Q = M                         # density rho = Q/M = 1

million = 1_000_000
# For a lattice containing M sites, define 1 sweep = M updates
num_sweeps = [1, 10, 100, 500, 1000, 2000, 5000]
snapshot_sequence = M * np.array(num_sweeps, dtype=np.int64)
snapshot_iterations = np.array(
    snapshot_sequence,
    dtype=np.int64,
)

boltz_array0 = initialize_boltzmann(Q, Lx, Ly)
assert boltz_array0.shape == (Lx, Ly)
assert boltz_array0.size == M

initial_total = boltz_array0.sum()

# Warm up Numba with a tiny 2-D run before timing the primary simulation.
warmup = np.ones((2, 5), dtype=np.int64)
warmup_checkpoints = np.array([10], dtype=np.int64)
_ = run_simulation_with_snapshots(warmup, warmup_checkpoints)

# Primary simulation loop. Keep at 1 for the snapshot study.
# Increase later to compare independent trajectories.
sim_runs = 1
run_times = []

for run in range(sim_runs):
    start_time = perf_counter()

    # Every run starts from the same pristine initial microstate.
    boltz_array = boltz_array0.copy()
    run_seed = boltz_seed + run

    print()
    print("CHECK: Beginning run:", run)
    print("CHECK: Initial top-left 2 x 5 block:")
    print(boltz_array[:2, :5])
    print("CHECK: Initial shape:", boltz_array.shape)
    print("CHECK: Initial total Q:", boltz_array.sum())

    boltz_array, evolved_snapshots = run_simulation_with_snapshots(
        boltz_array,
        snapshot_iterations,
        seed=run_seed,
    )

    # Include the exact initial condition as snapshot zero.
    snapshots = np.concatenate(
        (boltz_array0[np.newaxis, :, :], evolved_snapshots),
        axis=0,
    )
    # snapshot_labels = ["Initial"] + [f"{n:,} times" for n in snapshot_iterations]
    snapshot_labels = ["Initial"] + [f"{n:,} sweeps" for n in num_sweeps]

    print()
    print("CHECK: Final top-left 2 x 5 block:")
    print(boltz_array[:2, :5])
    print("CHECK: Final shape:", boltz_array.shape)
    print("CHECK: Final total Q:", boltz_array.sum())
    print()

    # Geometry and conservation checks
    assert boltz_array.shape == (Lx, Ly)
    assert boltz_array.size == M
    assert boltz_array.sum() == Q
    assert boltz_array.min() >= 0
    assert np.all(snapshots.sum(axis=(1, 2)) == Q)
    assert snapshots.min() >= 0

    elapsed_time = perf_counter() - start_time
    run_times.append(elapsed_time)

    print(f"Runtime: {elapsed_time:.3f} seconds")
    print(f"Initial total: {initial_total}")
    print(f"Final total:   {boltz_array.sum()}")
    print(f"Minimum value: {boltz_array.min()}")
    print(f"Maximum value: {boltz_array.max()}")

    # Plot all spatial snapshots with one common color scale.
    vmax = snapshots.max()
    fig, axes = plt.subplots(1, snapshots.shape[0], figsize=(18, 4), constrained_layout=True)

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
    fig.suptitle(f"One 2-D relaxation trajectory: {Lx} x {Ly}, Q = {Q:,}")


    # Occupation-number distribution of the final 2-D microstate.
    counts = np.bincount(boltz_array.ravel())
    values = np.arange(counts.size)

    # assume Poisson statistic on counts
    sigma_poisson = np.sqrt(counts)

    print("counts =", counts)
    print("values =", values)
    print("sigma_poisson =", np.sqrt(counts))
    print()

    # Strong conservation of quanta from the occupation distribution.
    assert (counts * values).sum() == Q
    assert counts.sum() == M

    occupied = counts > 0
    values_to_plot = values[occupied]
    counts_to_plot = counts[occupied]
    sigma_poisson_to_plot = sigma_poisson[occupied]

    fig2, ax = plt.subplots()
    ax.errorbar(
        values_to_plot,
        counts_to_plot,
        yerr=sigma_poisson_to_plot,
        marker="o",
        linestyle="None",
        capsize=3,
    )
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Number (n) of particles at a site")
    ax.set_ylabel("Number (N) of sites with n particles")
    ax.set_title("Final particle-count distribution")

plt.show()
