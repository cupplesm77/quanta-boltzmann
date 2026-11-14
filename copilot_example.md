## Current Model Example

As a current example the simulation:

- Lattice: 25 x 25
- Number of sites: M = 625
- Total quanta: Q = 625
- Mean occupancy: rho = Q/M = 1

The simulation begins from a deliberately low-entropy state: every site contains exactly one quantum. Random exchanges then drive the system toward statistical equilibrium.

An example snapshot sequence at

- 100            Get started
- M              1 Sweep
- 10 x M         10 Sweeps
- 100 M          100 Sweeps
- 1_000 x M      1000 Sweeps
- 1_500 x M      1500 Sweeps
- 2_000 x M      2000 Sweeps

iterations clearly shows the initial relaxation. By eye, however, the 100,000- and 1,000,000-iteration configurations are difficult to distinguish statistically.

This motivates a quantitative definition of equilibration.

---

## 1. Measure Time in Sweeps

For a lattice containing M sites, define

    1 sweep = M attempted updates

For the present 25 x 25 lattice,

    1 sweep = 625 iterations

Thus as an example, snapshots correspond approximately to:

| Iterations | Sweeps |
|---:|---:|
| 100 | 0.16 |
| 625 | 1.0 |
| 6250 | 10 |
| 62500 | 100 |
| 625000 | 1,000 |

Using sweeps rather than raw iteration counts will make comparisons between different lattice sizes more meaningful.

---

## 2. Equilibrium Is Statistical, Not Static

An equilibrium configuration does **not** stop changing.

Individual lattice occupancies continue to fluctuate. What becomes stationary are the statistical properties of the system.

Therefore:

> **Equilibrium means stationary statistics, not a stationary picture.**

Visual comparison of lattice snapshots is useful for intuition but is not a sufficient test for equilibrium.

---

## 3. Track the Occupancy Distribution

Define

    p_n(t) = N_n(t) / M

where N_n(t) is the number of sites containing n quanta at time t.

For the present system,

    rho = 1

and the large-system equilibrium prediction is the geometric distribution

    p_n(eq) = [1 / (1 + rho)] * [rho / (1 + rho)]^n

For rho = 1, this reduces to

    p_n(eq) = 1 / 2^(n+1)

Future development should compare the simulated p_n(t) with this theoretical prediction.

---

## 4. Define a Distance From Equilibrium

A useful single-number measure is the total variation distance

    D(t) = (1/2) * sum_n |p_n(t) - p_n(eq)|

During relaxation, D(t) should generally decrease.

At equilibrium it will not remain exactly zero because the finite lattice continues to fluctuate. Instead, it should reach a statistically stationary fluctuation level.

Plotting D(t) versus time, preferably on a logarithmic time scale, should provide a quantitative picture of equilibration.
