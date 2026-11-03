"""Two-Gaussian fit using KL(model || reference) for signed MCS angles.

Based on the objective and unit-variance constraint in Fruhwirth & Regler,
NIM A 456 (2001) 369-389, doi:10.1016/S0168-9002(00)00589-1, Sec. 4.
This is a regularized, BINNED approximation for Monte Carlo samples, not an
exact reproduction of their continuous-density calculation or optimizer.

Dependencies: NumPy and SciPy. Angles may be radians, mrad, or degrees;
returned sigmas have the same units as the input. Both means are fixed at 0.

Example:
    fit = fit_mixture_kl(theta_x, bins=200, pseudocount=0.5)
    weight = fit['core_weight']
    sigma_core = fit['sigma_core']
    sigma_tail = fit['sigma_tail']

For publication, check bins, coverage, pseudocount and sample-size dependence.
Bootstrap whole independent proton histories for sampling uncertainties;
the optimizer's inverse Hessian is not a sampling covariance here.
"""

import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr, rel_entr


def _normal_bin_probabilities(edges, variance):
    """Integrate N(0, variance) in each bin, retaining upper-tail precision."""
    lo = edges[:-1] / np.sqrt(variance)
    hi = edges[1:] / np.sqrt(variance)
    probability = np.where(
        lo >= 0,
        ndtr(-lo) - ndtr(-hi),  # Survival functions avoid 1 - 1 cancellation.
        ndtr(hi) - ndtr(lo),
    )
    return np.maximum(probability, 0.0)


def _fit_standardized_probabilities(edges, reference):
    """Fit normalized bin probabilities with total physical variance fixed at 1."""
    def model(parameters):
        epsilon, v_core = parameters
        v_tail = (1.0 - (1.0 - epsilon) * v_core) / epsilon
        return ((1.0 - epsilon) * _normal_bin_probabilities(edges, v_core) + epsilon * _normal_bin_probabilities(edges, v_tail))

    def objective(parameters):
        # The argument order implements the PAPER'S KL direction.
        return float(np.sum(rel_entr(model(parameters), reference)))

    bounds = [(1e-4, 0.4999), (1e-4, 0.9999)]  # epsilon, v_core
    starts = [(e, v) for e in (0.02, 0.08, 0.3) for v in (0.9, 0.6, 0.2)]
    results = [minimize(objective, start, method="L-BFGS-B", bounds=bounds,
            options={"maxiter": 2000, "maxls": 50, "ftol": 1e-12, "gtol": 1e-8}) 
                for start in starts]
    finite = [result for result in results if np.isfinite(result.fun)]
    if not finite:
        raise RuntimeError("All KL optimization attempts returned nonfinite costs.")
    best = min(finite, key=lambda result: result.fun)
    if not best.success:
        raise RuntimeError(f"Best KL optimization did not converge: {best.message}")
    at_bound = any(
        abs(value - low) < 1e-5 or abs(value - high) < 1e-5
        for value, (low, high) in zip(best.x, bounds))
    return best, model(best.x), at_bound


def fit_complex_kl(
    angles, bins=200, pseudocount=0.5, coverage=0.999, symmetrize=True,
):
    """Fit a centered core/tail mixture to signed projected scattering angles.

    Parameters
    ----------
    angles : 1-D array
        One signed projection per independent particle; finite values only.
        Do not supply unsigned polar angles. No recentering or outlier clipping
        is done. Remove a known incident-beam offset before calling if needed.
    bins : int
        Number of uniform interior bins; two additional overflow bins cover
        all angles outside the interior interval.
    pseudocount : float
        Positive count added to each bin to keep reference probabilities > 0.
        Default 0.5 is a regularization choice, NOT a value from the paper.
        Repeat with e.g. 0.1 and 1.0 to assess sensitivity.
    coverage : float, 0 < coverage <= 1
        Interior edges run from -q to q, where q is this quantile of abs(angle).
        All observations outside are retained in overflow bins. Large coverage
        can reduce core resolution at fixed bins. Changing this is a binning
        sensitivity check, not an angular cut.
    symmetrize : bool
        Pool opposite bins under the paper's zero-mean symmetric assumption.
        This does not double the sample size or recenter the beam.

    Returns
    -------
    dict
        core_weight (= 1 - epsilon), tail_weight (= epsilon), sigma_core,
        sigma_tail, v_core, v_tail, variance, kl and fit diagnostics.
        variance is mean(angle**2), the second moment about the assumed zero;
        it is the variance only when the population mean really is zero.
        Model bin probabilities include the overflow bins and sum to 1.
        The fitted full Gaussian mixture has total variance equal to variance.
    """
    theta = np.asarray(angles, dtype=float)
    if theta.ndim != 1 or theta.size < 20 or not np.all(np.isfinite(theta)):
        raise ValueError("angles must be a finite 1-D array with at least 20 entries.")
    if not (np.any(theta < 0) and np.any(theta > 0)):
        raise ValueError("Expected signed projected angles on both sides of zero.")
    if isinstance(bins, bool) or not isinstance(bins, (int, np.integer)) or bins < 10:
        raise ValueError("bins must be an integer >= 10.")
    if not np.isfinite(pseudocount) or pseudocount <= 0:
        raise ValueError("pseudocount must be finite and > 0.")
    if not np.isfinite(coverage) or not 0 < coverage <= 1:
        raise ValueError("coverage must be in (0, 1].")

    variance = float(np.mean(theta**2))
    if not np.isfinite(variance) or variance <= 0:
        raise ValueError("The second moment about zero must be finite and positive.")
    scale = np.sqrt(variance)
    x = theta / scale
    limit = float(np.quantile(np.abs(x), coverage))
    if not np.isfinite(limit) or limit <= 0:
        raise ValueError("The selected interior angular span is zero or nonfinite.")

    # Infinite end bins prevent silently dropping observations or model mass.
    edges = np.r_[-np.inf, np.linspace(-limit, limit, bins + 1), np.inf]
    counts = np.histogram(x, bins=edges)[0].astype(float)
    pooled_counts = 0.5 * (counts + counts[::-1]) if symmetrize else counts
    reference = pooled_counts + pseudocount
    reference /= reference.sum()

    best, model_probability, at_bound = _fit_standardized_probabilities(edges, reference)
    epsilon, v_core = map(float, best.x)
    v_tail = (1.0 - (1.0 - epsilon) * v_core) / epsilon
    pseudo_fraction = pseudocount * counts.size / (theta.size + pseudocount * counts.size)
    if pseudo_fraction > 0.01:
        warnings.warn("Pseudocounts supply over 1% of reference probability; "
                      "increase statistics or reduce the number of bins.", stacklevel=2)
    if abs(float(theta.mean())) > 5.0 * scale / np.sqrt(theta.size):
        warnings.warn("The sample mean is inconsistent with zero at about 5 SE; "
                      "check the incident reference direction.", stacklevel=2)
    if at_bound:
        warnings.warn("A mixture parameter is at an optimization bound; "
                      "inspect identifiability and the angular distribution.", stacklevel=2)
    if v_tail / v_core < 1.1:
        warnings.warn("The two widths nearly coincide; the tail weight is poorly "
                      "identified even if the fitted distribution is accurate.", stacklevel=2)

    return {
        "core_weight": 1.0 - epsilon,
        "tail_weight": epsilon,
        "sigma_core": float(scale * np.sqrt(v_core)),
        "sigma_tail": float(scale * np.sqrt(v_tail)),
        "v_core": v_core,
        "v_tail": v_tail,
        "variance": variance,
        "kl": float(best.fun),
        "success": bool(best.success),
        "at_bound": at_bound,
        "sample_mean": float(theta.mean()),
        "n_samples": int(theta.size),
        "pseudocount_fraction": float(pseudo_fraction),
        "overflow_fraction": float((counts[0] + counts[-1]) / theta.size),
        "edges": edges * scale,
        "counts": counts,
        "reference_probability": reference,
        "model_probability": model_probability,
    }


def fit_simple_kl(
    angles, bins=200, pseudocount=0.5, coverage=0.999, symmetrize=True,
    initial=None,
):
    """Return (core_weight, sigma_core, sigma_tail) for signed projected angles.

    Means are fixed at zero; sigmas retain the input units. The total second
    moment is fixed to mean(angles**2). Coverage sets the interior histogram
    span; overflow bins retain all remaining angles. Symmetrization pools
    opposite bins. Pseudocounts regularize empty bins: check their influence
    along with binning before interpreting tail weights.

    initial=(core_weight, sigma_core, sigma_tail) uses the input units.
    Its widths are rescaled together to the fixed second moment, preserving
    their ratio. A supplied guess uses one start; None uses nine default starts.
    """
    theta = np.asarray(angles, dtype=float)
    if theta.ndim != 1 or theta.size < 20 or not np.all(np.isfinite(theta)):
        raise ValueError("Provide at least 20 finite signed angles in a 1-D array.")
    if not (np.any(theta < 0) and np.any(theta > 0)):
        raise ValueError("Expected signed projected angles on both sides of zero.")
    if isinstance(bins, bool) or not isinstance(bins, (int, np.integer)) or bins < 10:
        raise ValueError("bins must be an integer >= 10.")
    if not np.isfinite(pseudocount) or pseudocount <= 0 or not 0 < coverage <= 1:
        raise ValueError("Require finite pseudocount > 0 and 0 < coverage <= 1.")

    scale = np.sqrt(np.mean(theta**2))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("The angular RMS must be finite and positive.")
    x = theta / scale
    limit = np.quantile(np.abs(x), coverage)
    if limit <= 0:
        raise ValueError("The selected histogram span is zero.")

    edges = np.r_[-np.inf, np.linspace(-limit, limit, bins + 1), np.inf]
    counts = np.histogram(x, bins=edges)[0].astype(float)
    if symmetrize:
        counts = 0.5 * (counts + counts[::-1])
    reference = counts + pseudocount
    reference /= reference.sum()

    def gaussian_bins(variance):
        lo, hi = edges[:-1] / np.sqrt(variance), edges[1:] / np.sqrt(variance)
        # Survival functions preserve precision in the positive tail.
        p = np.where(lo >= 0, ndtr(-lo) - ndtr(-hi), ndtr(hi) - ndtr(lo))
        return np.maximum(p, 0.0)

    def objective(parameters):
        epsilon, v_core = parameters
        v_tail = (1.0 - (1.0 - epsilon) * v_core) / epsilon
        model = (1.0 - epsilon) * gaussian_bins(v_core) + epsilon * gaussian_bins(v_tail)
        return float(rel_entr(model, reference).sum())

    if initial is None:
        starts = [(e, v) for e in (0.02, 0.08, 0.3) for v in (0.9, 0.6, 0.2)]
    else:
        guess = np.asarray(initial, dtype=float)
        if guess.shape != (3,) or not np.all(np.isfinite(guess)):
            raise ValueError("initial must be three finite values: (core_weight, sigma_core, sigma_tail).")
        w0, sc0, st0 = guess  # w0 is now the TAIL probability
        
        if not 0 < w0 < 0.5 or not 0 < sc0 < st0:
            raise ValueError("Require 0 < tail_weight < 0.5 and 0 < sigma_core < sigma_tail.")
        
        ratio_squared = (sc0 / st0)**2
        v0 = ratio_squared / ((1.0 - w0) * ratio_squared + w0)
        
        if not 1e-4 <= w0 <= 0.4999 or not 1e-4 <= v0 <= 0.9999:
            raise ValueError("The normalized initial guess is outside the optimizer bounds.")
        
        starts = [(w0, v0)]

    results = [
        minimize(
            objective, (e, v), method="L-BFGS-B",
            bounds=[(1e-4, 0.4999), (1e-4, 0.9999)],
            options={"maxiter": 2000, "maxls": 50, "ftol": 1e-12, "gtol": 1e-8},
        )
        for e, v in starts
    ]
    best = min(results, key=lambda r: r.fun if np.isfinite(r.fun) else np.inf)
    if not best.success or not np.isfinite(best.fun):
        raise RuntimeError(f"KL fit failed: {best.message}")

    epsilon, v_core = best.x
    v_tail = (1.0 - (1.0 - epsilon) * v_core) / epsilon
    return (
        float(epsilon),
        float(scale * np.sqrt(v_core)),
        float(scale * np.sqrt(v_tail)),
    )