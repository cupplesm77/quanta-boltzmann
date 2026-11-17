"""
autocorrelation_analysis.py

Autocorrelation and stationarity analysis for equilibrium-sampling data
produced by the statistical mechanics lattice simulation.

Purpose
-------
This script analyzes a CSV file containing measurements collected after
the simulation's designated equilibration period.

The analysis is intended to answer questions such as:

    1. Are successive recorded samples strongly correlated?

    2. How long is the characteristic autocorrelation time?

    3. How many effectively independent samples are contained in
       the recorded data?

    4. Do quantities such as D(t) and the population counts N_n appear
       stationary over the recorded equilibrium-sampling interval?

The script intentionally DOES NOT make a final binary declaration that
the system "has reached equilibrium." Instead, it reports quantitative
diagnostics that can be used to make that scientific judgment.

Typical usage
-------------
From the project root:

    python tools_analysis/autocorrelation_analysis.py equilibrium_samples.csv

The CSV file is expected to contain a sweep column and quantities such as:

    sweep
    D_t
    N_0
    N_1
    N_2
    ...

Author:
    Statistical Mechanics Simulation Project
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================================
# Configuration
# ============================================================================

# Maximum autocorrelation lag to display and analyze.
#
# This value is expressed in RECORDED SAMPLES, not simulation sweeps.
# For example, if samples are recorded every 5 sweeps, then lag 10
# corresponds to 50 simulation sweeps.
DEFAULT_MAX_LAG = 50

# Number of samples from the beginning and end of the recorded data that
# will be compared in the simple stationarity analysis.
#
# If the data set is too short for this many samples, the program will
# automatically use a smaller value.
DEFAULT_STATIONARITY_WINDOW = 50

# The D(t) column used by the current simulation output.
DEFAULT_D_COLUMN = "D_t"


# ============================================================================
# Data loading and validation
# ============================================================================

def load_equilibrium_data(csv_path: Path) -> pd.DataFrame:
    """
    Load equilibrium-sampling data from a CSV file.

    Parameters
    ----------
    csv_path : pathlib.Path
        Path to the CSV file produced by the simulation.

    Returns
    -------
    pandas.DataFrame
        Data frame containing the equilibrium measurements.

    Raises
    ------
    FileNotFoundError
        If the requested CSV file does not exist.

    ValueError
        If the file does not contain a 'sweep' column or contains too
        few rows for useful statistical analysis.
    """

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    data = pd.read_csv(csv_path)

    if "sweep" not in data.columns:
        raise ValueError(
            "The CSV file must contain a column named 'sweep'."
        )

    if len(data) < 3:
        raise ValueError(
            "At least three recorded samples are required for analysis."
        )

    return data


def determine_sample_interval(data: pd.DataFrame) -> float:
    """
    Determine the equilibrium sampling interval from the sweep column.

    The simulation normally records equilibrium samples at a constant
    interval. Rather than requiring the user to supply that interval
    manually, this function calculates it directly from the CSV data.

    Parameters
    ----------
    data : pandas.DataFrame
        Equilibrium measurements containing a 'sweep' column.

    Returns
    -------
    float
        Median number of simulation sweeps between recorded samples.

    Raises
    ------
    ValueError
        If the sweep values are not strictly increasing.
    """

    sweeps = data["sweep"].to_numpy(dtype=float)
    differences = np.diff(sweeps)

    if np.any(differences <= 0):
        raise ValueError(
            "Sweep values must be strictly increasing."
        )

    # Median is used because it is more robust than the mean if one
    # unusual spacing somehow appears in the CSV file.
    return float(np.median(differences))


def find_analysis_columns(
    data: pd.DataFrame,
    d_column: str = DEFAULT_D_COLUMN,
) -> list[str]:
    """
    Identify the simulation quantities to include in the analysis.

    The analysis includes D(t), when available, followed by all population
    columns whose names begin with 'N_'.

    Population columns are sorted numerically so that N_2 appears before
    N_10.

    Parameters
    ----------
    data : pandas.DataFrame
        Equilibrium measurements.

    d_column : str, optional
        Name of the column containing D(t).

    Returns
    -------
    list of str
        Column names selected for autocorrelation analysis.
    """

    columns = []

    if d_column in data.columns:
        columns.append(d_column)

    population_columns = [
        column for column in data.columns if column.startswith("N_")
    ]

    def population_number(column: str) -> int:
        """Return the integer portion of a population-column name."""
        try:
            return int(column.split("_", maxsplit=1)[1])
        except (IndexError, ValueError):
            return 10**9

    population_columns.sort(key=population_number)

    columns.extend(population_columns)

    return columns


# ============================================================================
# Autocorrelation calculations
# ============================================================================

def normalized_autocorrelation(
    values: np.ndarray,
    max_lag: int,
) -> np.ndarray:
    """
    Calculate the normalized autocorrelation function of a time series.

    For a measured quantity x, the normalized autocorrelation at lag k is

        rho(k) = C(k) / C(0)

    where C(k) is the covariance of measurements separated by k recorded
    samples.

    By definition:

        rho(0) = 1

    Parameters
    ----------
    values : numpy.ndarray
        One-dimensional sequence of recorded measurements.

    max_lag : int
        Largest lag, measured in recorded samples, to calculate.

    Returns
    -------
    numpy.ndarray
        Normalized autocorrelation values from lag 0 through max_lag.

    Notes
    -----
    The mean is calculated from the complete recorded sequence.

    The covariance at each lag is estimated from the overlapping sample
    pairs available at that lag.
    """

    x = np.asarray(values, dtype=float)

    if x.ndim != 1:
        raise ValueError("Autocorrelation input must be one-dimensional.")

    n = len(x)

    if n < 2:
        raise ValueError(
            "At least two values are required for autocorrelation."
        )

    max_lag = min(max_lag, n - 1)

    mean = np.mean(x)
    centered = x - mean

    variance = np.mean(centered * centered)

    # A constant series has zero variance, so its normalized
    # autocorrelation is mathematically undefined.
    if np.isclose(variance, 0.0):
        return np.full(max_lag + 1, np.nan)

    autocorrelation = np.empty(max_lag + 1, dtype=float)
    autocorrelation[0] = 1.0

    for lag in range(1, max_lag + 1):

        covariance = np.mean(
            centered[:-lag] * centered[lag:]
        )

        autocorrelation[lag] = covariance / variance

    return autocorrelation


def integrated_autocorrelation_time(
    autocorrelation: np.ndarray,
) -> float:
    """
    Estimate the integrated autocorrelation time.

    The estimator used here is

        tau_int = 1/2 + sum rho(k)

    where the summation proceeds through the initial sequence of positive
    autocorrelation values.

    Once rho(k) becomes zero or negative, the summation stops.

    This simple initial-positive-sequence rule prevents noisy,
    long-lag estimates from dominating the result.

    Parameters
    ----------
    autocorrelation : numpy.ndarray
        Normalized autocorrelation function beginning at lag zero.

    Returns
    -------
    float
        Estimated integrated autocorrelation time measured in recorded
        samples.

    Notes
    -----
    For essentially independent measurements, tau_int approaches 0.5.
    """

    if len(autocorrelation) == 0:
        return np.nan

    if np.isnan(autocorrelation[0]):
        return np.nan

    tau_int = 0.5

    for rho in autocorrelation[1:]:

        if not np.isfinite(rho):
            break

        if rho <= 0.0:
            break

        tau_int += rho

    return tau_int


def effective_sample_size(
    number_samples: int,
    tau_int: float,
) -> float:
    """
    Estimate the effective number of statistically independent samples.

    The standard relation is

        ESS = N / (2 * tau_int)

    where N is the number of recorded samples.

    Parameters
    ----------
    number_samples : int
        Number of recorded measurements.

    tau_int : float
        Integrated autocorrelation time measured in recorded samples.

    Returns
    -------
    float
        Estimated effective sample size.
    """

    if not np.isfinite(tau_int) or tau_int <= 0:
        return np.nan

    ess = number_samples / (2.0 * tau_int)

    # Numerical estimators can occasionally produce an ESS slightly above
    # the actual number of measurements. The effective number of samples
    # should not exceed the number actually recorded.
    return min(float(number_samples), ess)


# ============================================================================
# Stationarity analysis
# ============================================================================

def stationarity_statistics(
    values: np.ndarray,
    window: int,
) -> dict[str, float]:
    """
    Compare measurements near the beginning and end of the recorded run.

    This is intentionally a simple diagnostic rather than a formal
    equilibrium test. If the system is still relaxing appreciably,
    differences between the early and late means may reveal a systematic
    drift.

    Parameters
    ----------
    values : numpy.ndarray
        Recorded values of one simulation quantity.

    window : int
        Number of measurements to include in each comparison window.

    Returns
    -------
    dict
        Dictionary containing:

        - early_mean
        - late_mean
        - overall_mean
        - overall_std
        - difference

    Notes
    -----
    The early and late windows are chosen so that they do not overlap.
    """

    x = np.asarray(values, dtype=float)
    n = len(x)

    # Prevent the beginning and ending comparison windows from
    # overlapping.
    actual_window = min(window, n // 2)

    if actual_window < 1:
        raise ValueError(
            "Not enough measurements for stationarity analysis."
        )

    early = x[:actual_window]
    late = x[-actual_window:]

    early_mean = float(np.mean(early))
    late_mean = float(np.mean(late))
    overall_mean = float(np.mean(x))

    # ddof=1 gives the usual sample standard deviation when possible.
    overall_std = (
        float(np.std(x, ddof=1))
        if n > 1
        else 0.0
    )

    return {
        "early_mean": early_mean,
        "late_mean": late_mean,
        "overall_mean": overall_mean,
        "overall_std": overall_std,
        "difference": late_mean - early_mean,
    }


# ============================================================================
# Reporting
# ============================================================================

def analyze_quantity(
    data: pd.DataFrame,
    column: str,
    sample_interval_sweeps: float,
    max_lag: int,
    stationarity_window: int,
) -> dict:
    """
    Perform all analysis for one recorded simulation quantity.

    Parameters
    ----------
    data : pandas.DataFrame
        Complete equilibrium data set.

    column : str
        Name of the quantity to analyze.

    sample_interval_sweeps : float
        Number of simulation sweeps between recorded samples.

    max_lag : int
        Maximum autocorrelation lag to calculate.

    stationarity_window : int
        Number of early and late samples used for the stationarity
        comparison.

    Returns
    -------
    dict
        Dictionary containing autocorrelation, integrated autocorrelation
        time, effective sample size, and stationarity statistics.
    """

    values = data[column].to_numpy(dtype=float)

    acf = normalized_autocorrelation(
        values=values,
        max_lag=max_lag,
    )

    tau_samples = integrated_autocorrelation_time(acf)

    tau_sweeps = tau_samples * sample_interval_sweeps

    ess = effective_sample_size(
        number_samples=len(values),
        tau_int=tau_samples,
    )

    stationarity = stationarity_statistics(
        values=values,
        window=stationarity_window,
    )

    return {
        "column": column,
        "acf": acf,
        "rho_1": acf[1] if len(acf) > 1 else np.nan,
        "tau_samples": tau_samples,
        "tau_sweeps": tau_sweeps,
        "ess": ess,
        "stationarity": stationarity,
    }


def print_run_summary(
    data: pd.DataFrame,
    sample_interval_sweeps: float,
) -> None:
    """
    Print general information describing the recorded sampling run.
    """

    print()
    print("=" * 78)
    print("EQUILIBRIUM AUTOCORRELATION ANALYSIS")
    print("=" * 78)

    print(f"Recorded samples      : {len(data)}")
    print(f"Sample interval       : {sample_interval_sweeps:g} sweeps")
    print(f"First recorded sweep  : {data['sweep'].iloc[0]:g}")
    print(f"Last recorded sweep   : {data['sweep'].iloc[-1]:g}")
    print()


def print_autocorrelation_summary(results: list[dict]) -> None:
    """
    Print a compact table containing the principal autocorrelation results.
    """

    print("AUTOCORRELATION SUMMARY")
    print("-" * 78)

    header = (
        f"{'Quantity':<12}"
        f"{'rho(1)':>12}"
        f"{'tau(samples)':>16}"
        f"{'tau(sweeps)':>16}"
        f"{'ESS':>12}"
    )

    print(header)
    print("-" * 78)

    for result in results:
        print(
            f"{result['column']:<12}"
            f"{result['rho_1']:>12.4f}"
            f"{result['tau_samples']:>16.3f}"
            f"{result['tau_sweeps']:>16.3f}"
            f"{result['ess']:>12.1f}"
        )

    print()


def print_stationarity_summary(results: list[dict]) -> None:
    """
    Print early, late, and overall means for each analyzed quantity.
    """

    print("STATIONARITY SUMMARY")
    print("-" * 78)

    header = (
        f"{'Quantity':<12}"
        f"{'Early mean':>16}"
        f"{'Late mean':>16}"
        f"{'Overall mean':>16}"
        f"{'Late-Early':>16}"
    )

    print(header)
    print("-" * 78)

    for result in results:

        stats = result["stationarity"]

        print(
            f"{result['column']:<12}"
            f"{stats['early_mean']:>16.6g}"
            f"{stats['late_mean']:>16.6g}"
            f"{stats['overall_mean']:>16.6g}"
            f"{stats['difference']:>16.6g}"
        )

    print()


# ============================================================================
# Plotting
# ============================================================================

def plot_autocorrelations(
    results: list[dict],
    sample_interval_sweeps: float,
) -> None:
    """
    Plot the normalized autocorrelation functions.

    One plot is created for each analyzed quantity.

    The horizontal axis is expressed in simulation sweeps rather than
    recorded-sample number so that the plot corresponds directly to the
    simulation's time unit.

    Parameters
    ----------
    results : list of dict
        Analysis results returned by ``analyze_quantity``.

    sample_interval_sweeps : float
        Number of simulation sweeps between recorded samples.
    """

    for result in results:

        acf = result["acf"]

        lag_samples = np.arange(len(acf))
        lag_sweeps = lag_samples * sample_interval_sweeps

        fig, ax = plt.subplots(figsize=(8, 5))

        ax.plot(
            lag_sweeps,
            acf,
            marker="o",
            markersize=4,
            linewidth=1,
        )

        ax.axhline(
            0.0,
            color="black",
            linewidth=0.8,
        )

        ax.set_title(
            f"Autocorrelation: {result['column']}"
        )

        ax.set_xlabel("Lag (sweeps)")
        ax.set_ylabel("Normalized autocorrelation")

        ax.grid(True)

        fig.tight_layout()

    plt.show()


# ============================================================================
# Command-line interface
# ============================================================================

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Perform autocorrelation and stationarity analysis on "
            "equilibrium-sampling data from the StatMech simulation."
        )
    )

    parser.add_argument(
        "csv_file",
        type=Path,
        help="Path to the equilibrium CSV file.",
    )

    parser.add_argument(
        "--max-lag",
        type=int,
        default=DEFAULT_MAX_LAG,
        help=(
            "Maximum autocorrelation lag in recorded samples "
            f"(default: {DEFAULT_MAX_LAG})."
        ),
    )

    parser.add_argument(
        "--stationarity-window",
        type=int,
        default=DEFAULT_STATIONARITY_WINDOW,
        help=(
            "Number of early and late samples used for the "
            "stationarity comparison "
            f"(default: {DEFAULT_STATIONARITY_WINDOW})."
        ),
    )

    return parser.parse_args()


# ============================================================================
# Main
# ============================================================================

def main() -> None:
    """
    Run the equilibrium autocorrelation analysis.
    """

    args = parse_arguments()

    # ------------------------------------------------------------------------
    # Load the data.
    # ------------------------------------------------------------------------

    data = load_equilibrium_data(args.csv_file)

    # ------------------------------------------------------------------------
    # Determine the actual number of recorded measurements.
    # We obtain this from the data itself rather than assuming a value.
    # ------------------------------------------------------------------------

    sample_interval_sweeps = determine_sample_interval(data)

    # ------------------------------------------------------------------------
    # Find the quantities to analyze.
    # ------------------------------------------------------------------------

    analysis_columns = find_analysis_columns(data)

    if not analysis_columns:
        raise ValueError(
            "No D_t or N_n columns were found in the CSV file."
        )

    # ------------------------------------------------------------------------
    # Perform the analysis.
    # ------------------------------------------------------------------------

    results = []

    for column in analysis_columns:

        result = analyze_quantity(
            data=data,
            column=column,
            sample_interval_sweeps=sample_interval_sweeps,
            max_lag=args.max_lag,
            stationarity_window=args.stationarity_window,
        )

        results.append(result)

    # ------------------------------------------------------------------------
    # Print numerical summaries.
    # ------------------------------------------------------------------------

    print_run_summary(
        data=data,
        sample_interval_sweeps=sample_interval_sweeps,
    )

    print_autocorrelation_summary(results)

    print_stationarity_summary(results)

    # ------------------------------------------------------------------------
    # Plot the autocorrelation functions.
    # ------------------------------------------------------------------------

    plot_autocorrelations(
        results=results,
        sample_interval_sweeps=sample_interval_sweeps,
    )


if __name__ == "__main__":
    main()