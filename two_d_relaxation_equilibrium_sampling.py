"""Simulate relaxation and equilibrium fluctuations on a 2-D Einstein solid.

The lattice contains ``M = Lx * Ly`` sites and a conserved total of ``Q``
discrete energy quanta. A Monte Carlo update selects a donor site uniformly;
when the donor is occupied, one quantum is transferred to a uniformly selected
receiver. Donor and receiver may be anywhere on the lattice, so the model has
global (mean-field) rather than nearest-neighbor exchange.

The script performs two complementary experiments:

1. Independent relaxation runs of several durations, recording the
   total-variation distance D(t) after every sweep.
2. One long equilibrium run, recording occupation-bin counts at fixed sweep
   intervals for later fluctuation and autocorrelation analysis.

Raw equilibrium observations are written to ``equilibrium_bin_counts.csv``.
Statistical analysis is intentionally kept outside this simulation so these
results can be reanalyzed without repeating the expensive Monte Carlo run.

Notes
-----
A sweep is defined as ``M`` attempted single-quantum transfers. Random seeds
are fixed for reproducibility. The primary numerical kernels use Numba's
``@njit`` decorator and therefore incur compilation overhead on their first
call in a Python process.
"""

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
    """Calculate total-variation distance from geometric equilibrium.

    Parameters
    ----------
    boltz_array : numpy.ndarray
        Integer lattice whose entries are site occupation numbers. The array is
        read only and may have any shape.
    rho : float
        Mean occupation ``Q / M``. It determines the equilibrium geometric
        distribution ``p_n = (1 / (1 + rho)) * (rho / (1 + rho))**n``.

    Returns
    -------
    float
        Total-variation distance between the measured occupation-frequency
        distribution and the infinite-support geometric distribution. The
        exact theoretical tail above the largest observed occupation is
        included in the result.

    Notes
    -----
    The historical function name contains ``equillibrium``. It is retained to
    avoid changing existing calls and saved study material.
    """
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
    """Run one relaxation experiment and record D(t) after every sweep.

    Parameters
    ----------
    boltz_array : numpy.ndarray
        Initial integer lattice. It is modified in place and returned in its
        final state.
    total_sweeps : int
        Number of sweeps to simulate. One sweep equals one attempted update per
        lattice site, or ``boltz_array.size`` attempted transfers.
    seed : int
        Random-number seed for this independent trajectory.
    rho : float
        Conserved mean occupation ``Q / M`` used to construct the theoretical
        geometric equilibrium distribution.

    Returns
    -------
    boltz_array : numpy.ndarray
        Final lattice microstate after ``total_sweeps`` sweeps.
    distances : numpy.ndarray
        One-dimensional float array of length ``total_sweeps + 1``. Element 0
        is D(0); element ``s`` is the distance after sweep ``s``.

    Notes
    -----
    The distance calculation is written explicitly inside this Numba kernel.
    Calling the separate Python helper would prevent the complete time-critical
    loop from remaining straightforwardly JIT compiled.
    """
    # Seed inside the compiled routine so a run is reproducible in isolation.
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


@njit
def run_equilibrium_sampling(
        boltz_array,
        equilibration_sweeps,
        sample_interval_sweeps,
        num_samples,
        seed,
        rho,
        max_recorded_n,
):
    """Equilibrate one lattice and sample occupation counts at fixed intervals.

    The relaxation transient is not stored. After the chosen equilibration
    period, the same isolated lattice keeps evolving and is measured every
    ``sample_interval_sweeps`` sweeps. Each measurement records the absolute
    sweep number, D(t), counts N_0 through N_max, and an overflow count.

    Parameters
    ----------
    boltz_array : numpy.ndarray
        Initial integer lattice. It is modified in place.
    equilibration_sweeps : int
        Number of complete sweeps discarded before equilibrium sampling begins.
    sample_interval_sweeps : int
        Sweep spacing between successive saved measurements.
    num_samples : int
        Number of equilibrium measurements to collect.
    seed : int
        Random-number seed for the complete equilibration-and-sampling run.
    rho : float
        Mean occupation ``Q / M`` used for the theoretical equilibrium law.
    max_recorded_n : int
        Largest occupation with its own saved count column. Sites above this
        value are included in the separate overflow count.

    Returns
    -------
    boltz_array : numpy.ndarray
        Final lattice microstate.
    sweep_numbers : numpy.ndarray
        Integer sweep number associated with each saved measurement.
    distances : numpy.ndarray
        Total-variation distance D(t) at each saved measurement.
    occupation_history : numpy.ndarray
        Integer array of shape ``(num_samples, max_recorded_n + 1)``. Entry
        ``[i, n]`` is the number of sites with occupation ``n`` in sample ``i``.
    overflow_history : numpy.ndarray
        Number of sites above ``max_recorded_n`` in each sample.

    Notes
    -----
    D(t) is calculated from the complete observed occupation range, independent
    of ``max_recorded_n``. Thus the storage limit cannot truncate the distance
    calculation. Keeping an explicit overflow count also prevents silent loss
    of sites when an unexpectedly large occupation occurs.
    """
    np.random.seed(seed)
    flat = boltz_array.ravel()
    M = flat.size

    sweep_numbers = np.empty(num_samples, dtype=np.int64)
    distances = np.empty(num_samples, dtype=np.float64)
    occupation_history = np.zeros(
        (num_samples, max_recorded_n + 1), dtype=np.int64
    )
    overflow_history = np.empty(num_samples, dtype=np.int64)

    ratio = rho / (1.0 + rho)
    prefactor = 1.0 / (1.0 + rho)
    final_sweep = equilibration_sweeps + sample_interval_sweeps * num_samples
    sample_index = 0

    for sweep in range(1, final_sweep + 1):
        # One sweep = M attempted single-quantum transfers.
        for _ in range(M):
            index1 = np.random.randint(M)
            if flat[index1] > 0:
                index2 = np.random.randint(M)
                if index2 != index1:
                    flat[index1] -= 1
                    flat[index2] += 1

        # Do not measure the relaxation transient.  Sampling begins only after
        # the chosen equilibration period and then proceeds at a fixed spacing.
        if sweep > equilibration_sweeps and \
                (sweep - equilibration_sweeps) % sample_interval_sweeps == 0:
            counts = np.zeros(max_recorded_n + 1, dtype=np.int64)
            overflow = 0
            for value in flat:
                if value <= max_recorded_n:
                    counts[value] += 1
                else:
                    overflow += 1

            # Compute D(t) exactly from the same equilibrium microstate.
            # This calculation is intentionally independent of max_recorded_n,
            # which is only a storage limit for the CSV occupation columns.
            sample_max_n = flat.max()
            distance_counts = np.zeros(sample_max_n + 1, dtype=np.int64)
            for value in flat:
                distance_counts[value] += 1

            distance = 0.0
            for n in range(distance_counts.size):
                p_measured = distance_counts[n] / M
                p_equilibrium = prefactor * ratio**n
                distance += abs(p_measured - p_equilibrium)
            distance += ratio**distance_counts.size

            sweep_numbers[sample_index] = sweep
            distances[sample_index] = 0.5 * distance
            occupation_history[sample_index, :] = counts
            overflow_history[sample_index] = overflow
            sample_index += 1

    return (
        boltz_array,
        sweep_numbers,
        distances,
        occupation_history,
        overflow_history,
    )

# =============================================================================
# User-adjustable experiment configuration
# =============================================================================
# Lattice geometry and conserved total. Large defaults are intentional: they
# suppress finite-size noise but require substantial memory and runtime.
Lx = 1000
Ly = 1000
M = Lx * Ly
Q = M                         # density rho = Q/M = 1
rho = Q / M

# Each value is now an independent relaxation experiment.
# num_sweeps = [1, 10, 100, 500, 1000, 2000, 5000]
num_sweeps = [1, 10, 20, 30, 40, 50, 100]

# Construct one pristine reference state. Every independent relaxation case
# receives a copy so case duration, rather than inherited history, is compared.
boltz_array0 = initialize_boltzmann(Q, Lx, Ly)
assert boltz_array0.shape == (Lx, Ly)
assert boltz_array0.size == M
initial_total = boltz_array0.sum()

# Warm up the JIT-compiled relaxation routine before timing experiments.
warmup = np.ones((2, 5), dtype=np.int64)
_ = run_relaxation_history(warmup, 1, boltz_seed, 1.0)

# Results remain in memory for plotting after all independent runs complete.
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


# -----------------------------------------------------------------------------
# Equilibrium-fluctuation experiment
# -----------------------------------------------------------------------------
# Keep this experiment independent of the relaxation/snapshot cases above.
# Begin again from the same pristine low-entropy state, allow the lattice to
# equilibrate, and only then start collecting N_n measurements.
#
# These are intentionally plain, user-adjustable controls.  In particular,
# sample_interval_sweeps can later be changed after we study autocorrelation.
equilibration_sweeps = num_sweeps[-1]
sample_interval_sweeps = 10
num_equilibrium_samples = 500
max_recorded_n = 64
sampling_seed = boltz_seed + 1000

equilibrium_array = boltz_array0.copy()
print()
print("CHECK: Beginning equilibrium-fluctuation sampling experiment")
print(f"Equilibration period: {equilibration_sweeps:,} sweeps")
print(f"Sampling interval:    {sample_interval_sweeps:,} sweep(s)")
print(f"Samples requested:    {num_equilibrium_samples:,}")

sampling_start = perf_counter()
(
    equilibrium_array,
    equilibrium_sample_sweeps,
    equilibrium_sample_distances,
    equilibrium_counts,
    equilibrium_overflow,
) = run_equilibrium_sampling(
    equilibrium_array,
    equilibration_sweeps,
    sample_interval_sweeps,
    num_equilibrium_samples,
    sampling_seed,
    rho,
    max_recorded_n,
)
sampling_elapsed = perf_counter() - sampling_start

# Conservation checks are kept close to data generation so a corrupted run is
# rejected before its measurements are written to disk.
assert equilibrium_array.sum() == Q
assert equilibrium_array.min() >= 0
assert np.all(equilibrium_counts.sum(axis=1) + equilibrium_overflow == M)

# Write RAW equilibrium observations only.  Statistical comparisons (mean,
# variance, Fano factor, autocorrelation, Poisson tests) belong in a separate
# analysis program so the simulation remains a clean data-producing engine.
#
# Trim unused high-n columns while retaining every occupation actually observed
# when it lies within max_recorded_n.  The overflow column remains as a guard.
occupied_columns = np.where(equilibrium_counts.sum(axis=0) > 0)[0]
highest_recorded_n = int(occupied_columns[-1]) if occupied_columns.size else 0
counts_to_save = equilibrium_counts[:, :highest_recorded_n + 1]

output_data = np.column_stack((
    equilibrium_sample_sweeps,
    equilibrium_sample_distances,
    counts_to_save,
    equilibrium_overflow,
))
header_columns = ["sweep", "D_t"]
header_columns += [f"N_{n}" for n in range(highest_recorded_n + 1)]
header_columns += [f"N_gt_{max_recorded_n}"]
np.savetxt(
    "equilibrium_bin_counts.csv",
    output_data,
    delimiter=",",
    header=",".join(header_columns),
    comments="",
    fmt=["%d", "%.12g"] + ["%d"] * (counts_to_save.shape[1] + 1),
)

print(f"Sampling runtime:     {sampling_elapsed:.3f} seconds")
print(f"Samples written:      {num_equilibrium_samples:,}")
print(f"Highest saved n:      {highest_recorded_n}")
print(f"Maximum overflow:     {equilibrium_overflow.max():,} site(s)")
print("Data file:            equilibrium_bin_counts.csv")

# =============================================================================
# Diagnostic graphics
# =============================================================================
# These plots are exploratory study aids rather than publication-formatted
# figures. They display spatial relaxation, D(t), and the final count spectrum.
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
# Use one shared color scale so color has the same meaning in every panel.
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
# sqrt(N) error bars provide a familiar Poisson counting-noise reference.
# They are not, by themselves, a proof that equilibrium bin fluctuations are
# Poisson; that question belongs to the separate time-series analysis.
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
