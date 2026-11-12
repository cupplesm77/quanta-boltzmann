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
    """Run one trajectory and save the 2-D state at selected iterations."""
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
    """Distribute Q quanta approximately uniformly over an Lx by Ly array."""
    M = Lx * Ly
    q, r = divmod(Q, M)

    # Give every lattice point q quanta.
    boltzmann = np.full((Lx, Ly), q, dtype=np.int64)

    # Spread any remaining quanta approximately uniformly using a flat view.
    if r > 0:
        flat = boltzmann.ravel()
        indices = (np.arange(r) * M) // r
        flat[indices] += 1

    return boltzmann


# Lattice geometry and conserved total
Lx = 25
Ly = 25
M = Lx * Ly
Q = M                         # density rho = Q/M = 1

million = 1_000_000
snapshot_iterations = np.array(
                              [1_00, 0.001 * million, 0.01 * million, 0.1 * million, million],
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
    snapshot_labels = ["Initial"] + [f"{n:,} iterations" for n in snapshot_iterations]

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

    print("counts =", counts)
    print("values =", values)
    print()

    # Strong conservation of quanta from the occupation distribution.
    assert (counts * values).sum() == Q
    assert counts.sum() == M

    occupied = counts > 0
    values_to_plot = values[occupied]
    counts_to_plot = counts[occupied]

    fig2, ax = plt.subplots()
    ax.plot(values_to_plot, counts_to_plot, marker="o", linestyle="None")
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Number (n) of particles at a site")
    ax.set_ylabel("Number (N) of sites with n particles")
    ax.set_title("Final particle-count distribution")

plt.show()
