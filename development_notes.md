# StatMech Simulation: Equilibration Development Plan

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

---

## 5. Test Several Initial Conditions

Eventually repeat the experiment using very different starting configurations, for example:

1. Uniform distribution
2. Highly concentrated distribution
3. Random distribution
4. Spatially clustered distribution

If these different initial states all evolve toward the same equilibrium statistics, this provides strong evidence 
that the simulation has lost memory of its initial condition.

---

## 6. Later: Measure Autocorrelation

Two questions must eventually be distinguished:

1. **Equilibration time:** How long should the initial transient ("burn-in") be discarded?
2. **Sampling interval:** Once equilibrium has been reached, how widely separated should measurements be to obtain 
substantially independent samples?

The second question requires investigation of autocorrelation times.

---

## 7. Investigate Lattice-Size Scaling

The number of iterations required for equilibration should not be assumed to depend only on the number of lattice sites.

Future experiments should systematically vary lattice size and measure equilibration in sweeps.

This will allow the simulation itself to reveal how relaxation time scales with system size for the particular update rule being used.

---

## Immediate Next Development Steps

When simulation development resumes:

1. Compute p_n(t) at selected times.
2. Compare p_n(t) with the predicted geometric equilibrium distribution.
3. Calculate and plot D(t).
4. Estimate an equilibration/burn-in time from the stationary region of D(t).
5. Repeat from different initial conditions.
6. Investigate autocorrelation.
7. Repeat for different lattice sizes and study scaling.

---

## Current Stopping Point

The simulation has reached a useful stage:

- Low-to-high entropy relaxation is visually evident.
- The expected equilibrium occupancy distribution is known analytically.
- The limitations of judging equilibrium from snapshots alone are understood.
- A quantitative plan for determining equilibration has been established.

Simulation development can therefore pause here.

The next priority is to work systematically through the previously developed StatMech exercises, using them to strengthen intuition about multiplicity, probability, equilibrium, fluctuations, and the emergence of the geometric occupancy distribution.

The simulation can then be extended naturally as those concepts are developed.

ASCII convention going forward

For this StatMech work, I think this is a good convention for our working Markdown:

rho = Q/M
p_n
p_n(t)
p_n(eq)
p_n(eq) = [1/(1+rho)] * [rho/(1+rho)]^n
D(t) = (1/2) * sum_n |p_n(t) - p_n(eq)|
M = Lx * Ly


That should keep the notes immediately readable both as raw Markdown and when rendered, while still making the 
mathematical meaning unambiguous.

The substance is unchanged; the main change is replacing display LaTeX with compact ASCII notation suitable for 
your simulation development notes.
