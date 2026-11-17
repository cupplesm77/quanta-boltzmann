Statistical Mechanics Analysis Tools
This directory contains analysis tools used with the statistical mechanics lattice simulation.

The simulation code is responsible for generating physical data. The tools in this directory are intended to analyze 
that data separately, allowing the simulation and the statistical analysis to remain logically independent.

The first tool in this directory is:

autocorrelation_analysis.py
It analyzes equilibrium-sampling data produced by the lattice simulation.

1. Purpose
A simulation may appear visually or numerically close to equilibrium without that observation alone establishing 
2. that measurements taken from the system are statistically independent.

The autocorrelation analysis tool is intended to help answer several questions:

Does the recorded portion of the run appear stationary?
How strongly is each recorded quantity correlated with earlier values?
What is the integrated autocorrelation time?
Is the chosen sampling interval long enough to provide approximately independent measurements?
How many effectively independent measurements are represented by the recorded samples?
The program intentionally does not make a binary declaration such as:

EQUILIBRIUM ACHIEVED: YES
Instead, it provides quantitative diagnostics that can be used to make that scientific judgment.

2. Current Analysis
The program currently analyzes:

D_t
N_0
N_1
N_2
...
where:

D_t is the simulation's equilibrium-distance diagnostic.
N_n is the number of lattice sites having occupation number n.
The program automatically discovers the available N_n columns in the CSV file. The maximum recorded occupation 
number therefore does not need to be hard-coded into the analysis script.

3. Expected Input
The program expects a CSV file generated during the equilibrium-sampling portion of the simulation.

The CSV must contain a column named:

sweep
and normally contains additional columns such as:

sweep
D_t
N_0
N_1
N_2
N_3
...
The sweep column must increase monotonically.

The sampling interval is determined directly from the sweep values in the CSV rather than being supplied separately 
to the analysis program.

For example, sweep values:

75
80
85
90
...
indicate:

sample_interval_sweeps = 5
4. Running the Program
From the project root directory:

python analysis_tools/autocorrelation_analysis.py <csv_file>
Example:

python analysis_tools/autocorrelation_analysis.py equilibrium_bin_counts_70_5_300.csv
Optional command-line arguments are available.

For example:

python analysis_tools/autocorrelation_analysis.py equilibrium_bin_counts_70_5_300.csv --max-lag 50 --stationarity-window 50
The defaults are:

max_lag = 50 samples
stationarity_window = 50 samples
5. Output
The program produces two principal numerical summaries:

Autocorrelation summary
Stationarity summary
It also produces autocorrelation plots for selected important observables.

6. Autocorrelation
For a recorded quantity x, the normalized autocorrelation function is calculated as a function of lag.

Conceptually:

rho(k) = C(k) / C(0)
where:

C(k) = covariance between measurements separated by lag k
C(0) = variance of the quantity
Therefore:

rho(0) = 1
A value of rho(k) near zero indicates little correlation between measurements separated by that lag.

A large positive value indicates that measurements retain significant memory of earlier measurements.

7. Lag and Simulation Time
Autocorrelation lag is calculated in recorded samples.

The program also converts the result into the simulation time unit of sweeps.

For example, if:

sample_interval_sweeps = 5
then:

lag 1 = 5 sweeps
lag 2 = 10 sweeps
lag 3 = 15 sweeps
...
This distinction is important when interpreting autocorrelation results.

8. Integrated Autocorrelation Time
The program estimates the integrated autocorrelation time:

tau_int = 1/2 + sum rho(k)
using the initial sequence of positive autocorrelation values.

The summation stops when the autocorrelation becomes zero or negative.

The result is first calculated in recorded samples and then converted into simulation sweeps.

The output therefore contains both:

tau(samples)
tau(sweeps)
For nearly independent recorded measurements:

tau(samples) ~ 0.5
If samples are recorded every 5 sweeps, this corresponds to:

tau(sweeps) ~ 2.5
9. Effective Sample Size
Successive simulation measurements may not be fully independent.

The effective sample size (ESS) estimates how many independent measurements the recorded data approximately represent.

The program uses:

ESS = N / (2 * tau_int)
where:

N       = number of recorded samples
tau_int = integrated autocorrelation time in sample units
For nearly independent data:

ESS ~ N
Strong autocorrelation reduces the effective sample size.

10. Stationarity Analysis
The program also performs a simple stationarity diagnostic.

It compares a window near the beginning of the recorded data with a window near the end.

The report includes:

Early mean
Late mean
Overall mean
Late - Early
With the default settings:

stationarity_window = 50
the first 50 and last 50 samples are compared, assuming the data set is large enough to support non-overlapping windows.

A system that has reached a stationary regime should fluctuate around stable statistical values rather than show 
continued systematic drift.

This comparison is a diagnostic, not by itself a formal proof of equilibrium.

11. Autocorrelation Plots
Numerical autocorrelation statistics are calculated for all available recorded population bins.

For clarity, plots are currently limited to:

D_t
N_0
N_1
N_2
N_3
N_4
N_5
N_6
This avoids generating large numbers of plots for high-order population bins that may contain very few observations.

The plots display normalized autocorrelation versus lag measured in simulation sweeps.

12. Zero-Variance Quantities and NaN
Some high-occupation bins may contain zero particles during every recorded sample.

For example, a column might contain:

N_20 = 0, 0, 0, 0, ...
Such a sequence has:

variance = 0
Normalized autocorrelation requires division by the variance. Consequently, autocorrelation is mathematically 
undefined for a constant sequence.

The program reports these results as:

nan
This is expected behavior and does not indicate a program failure.

13. Caution With Sparsely Populated Bins
High-order occupation bins can contain very few events.

For example, a bin whose mean population is much less than one may consist primarily of zeros with occasional 
nonzero counts.

Autocorrelation estimates for such sparse quantities can fluctuate substantially and should be interpreted cautiously.

For equilibrium assessment, greater diagnostic weight should generally be placed on well-populated observables 
together with D_t, rather than on a single sparsely populated high-order bin.

14. Initial Validation Case
The program was initially tested using a 100 x 100 lattice run with:

equilibration_sweeps      = 70
sample_interval_sweeps    = 5
num_equilibrium_samples   = 300
max_recorded_n            = 26
The resulting CSV contained:

Recorded samples      : 300
Sample interval       : 5 sweeps
First recorded sweep  : 75
Last recorded sweep   : 1570
For this validation case, selected results were:

Quantity    rho(1)    tau(samples)    tau(sweeps)    ESS

D_t        -0.0127       0.500           2.500       300.0
N_0         0.0917       0.592           2.958       253.5
N_1         0.0120       0.512           2.560       293.0
N_2         0.0903       0.663           3.317       226.1
N_3         0.0397       0.540           2.699       277.9
N_4        -0.0200       0.500           2.500       300.0
N_5         0.0592       0.585           2.927       256.2
N_6         0.1022       0.615           3.077       243.8
The principal observables therefore showed relatively short autocorrelation times compared with the 5-sweep 
sampling interval.

The stationarity analysis also showed little early-to-late change in the principal observables.

For example:

D_t early mean = 0.00731861
D_t late mean  = 0.00732899
These results were interpreted as evidence that the 100 x 100 system was stationary when equilibrium sampling 
began following the 70-sweep equilibration period.

Important:

The equilibrium CSV begins at sweep 75. It therefore demonstrates the behavior of the system after the designated 
70-sweep equilibration period. It does not determine whether equilibrium may have been reached earlier than 70 sweeps.

15. Scientific Interpretation
The analysis separates two related but distinct questions:

Question 1:
Has the system stopped showing systematic relaxation?

Question 2:
How statistically independent are measurements taken from the
apparently equilibrated system?
The stationarity analysis primarily addresses the first question.

The autocorrelation analysis primarily addresses the second.

Neither diagnostic should be regarded in isolation as an absolute proof of thermodynamic equilibrium.

16. Design Philosophy
The analysis tools are intentionally maintained separately from the main simulation code.

The working philosophy is:

Simulation code:
    Generate the physical evolution and record data.

Analysis code:
    Examine the recorded data and quantify its statistical behavior.

Scientific interpretation:
    Use the diagnostics to determine what conclusions are justified.
This separation helps prevent analysis assumptions from becoming hidden inside the simulation itself.

It also permits analysis methods to evolve without unnecessarily modifying the underlying simulation.

17. Current Limitations
The current program is intentionally conservative.

It does not currently:

Declare equilibrium automatically.
Determine the minimum number of equilibration sweeps.
Select an optimal sampling interval automatically.
Perform statistical hypothesis tests for stationarity.
Determine confidence intervals for autocorrelation estimates.
Automatically classify sparse high-order bins.
Compare multiple simulation runs simultaneously.
These capabilities may be considered later as the behavior of the simulation is studied across additional 
lattice sizes and parameter choices.

18. Possible Future Development
Potential extensions include:

relaxation_analysis.py
equilibrium_analysis.py
distribution_analysis.py
Possible future capabilities include:

Automated comparison of several equilibration times.
Relaxation-time estimation.
Statistical comparison of early and late equilibrium windows.
Confidence intervals.
Automated sparse-bin handling.
Comparison of autocorrelation behavior with lattice size.
Automated CSV or text summary reports.
Saving analysis plots directly to an output directory.
These features should be added incrementally as the physical and statistical requirements become better established.

19. Dependencies
The current autocorrelation analysis uses:

Python
NumPy
pandas
Matplotlib
The required Python packages must be available in the active Python environment.

20. Project Principle
This analysis package follows a deliberately conservative development rule:

Measure first, automate the measurement second, and establish decision rules only after the statistical behavior 
of the simulation is well understood.

The objective is not merely to obtain an answer from the simulation, but to understand why the available numerical 
evidence justifies that answer.