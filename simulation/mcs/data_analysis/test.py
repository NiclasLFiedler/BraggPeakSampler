import operator
import numpy as np

def _beta(E0, T):
    return np.sqrt(T * (T + 2 * E0) / (T + E0)**2)

def _momentum(E0, T):
    return np.sqrt(T * (T + 2 * E0))

def _energy_at_depth(z, R0, E0, p_exp):
    return E0 * (1-z/R0)**(1/p_exp)

def two_grid(R0, switch_depth, coarse_step, fine_step): ##returns depths with 2 widths and dz
    """Return interior reach depths z and integration widths dz, starting at 0.

    Steps are maximum widths: each region uses equally spaced intervals no
    wider than requested. Includes the switch exactly and excludes R0.
    There is a final unresolved interval (z[-1], R0), at most fine_step wide.
    """
    values = np.asarray([R0, switch_depth, coarse_step, fine_step], float)
    
    if not np.all(np.isfinite(values)):
        raise ValueError('Grid parameters must be finite.')
    if not (0 < switch_depth < R0 and coarse_step > 0 and fine_step > 0):
        raise ValueError('Require 0 < switch_depth < R0 and positive steps.')
    
    n_front = int(np.ceil(switch_depth / coarse_step))
    n_tail = int(np.ceil((R0-switch_depth) / fine_step))
    edges = np.r_[np.linspace(0.0, switch_depth, n_front+1),
                  np.linspace(switch_depth, R0, n_tail+1)[1:]]
    return edges[1:-1], np.diff(edges)[:-1]


def _variance(V): ##checks variance
    V = np.asarray(V, dtype=float)
    if V.ndim != 1 or V.size == 0 or not np.all(np.isfinite(V)):
        raise ValueError('V must be a nonempty finite 1D array in rad^2.')
    if np.any(V < 0) or np.any(np.diff(V) < 0):
        raise ValueError('V must be nonnegative and nondecreasing. Do not '
                         'silently clip negative measured increments.')
    return V


def _grid(z, dz, R0):
    z, dz = np.asarray(z, float), np.asarray(dz, float)
    if (z.ndim != 1 or z.size == 0 or dz.shape != z.shape
            or not np.all(np.isfinite(z)) or not np.all(np.isfinite(dz))
            or not np.isfinite(R0) or R0 <= 0):
        raise ValueError('Invalid grid or R0.')
    if np.any(dz <= 0) or np.any(np.diff(z) <= 0) or np.any(z >= R0):
        raise ValueError('Need positive intervals and increasing depths < R0.')
    if not np.allclose(z, np.cumsum(dz), rtol=1e-12, atol=1e-13):
        raise ValueError('dz must cover every interval from 0 to each depth.')
    return z, dz


def simulation_grid(depths, cumulative_variance, R0, layer_ids=None):
    """Use actual downstream boundary depths and measured V without interpolation.

    depths: one depth per consecutive scored layer in physical order, cm.
    cumulative_variance: matching ONE-projection variances in rad^2.
    R0: fixed reference path range, cm.
    layer_ids: optional matching IDs, required to be 0,1,... . Providing these
        detects omitted layers (including missing front layers). Without IDs,
        completeness must be established by the caller against the geometry.

    Only depths strictly below R0 are used. Nonfinite variances at/beyond R0
    are irrelevant; nonfinite/decreasing variances below R0 raise. There must
    be no gaps of unmodelled material or omitted layers in the input prefix.
    A depth equal to 0 is not an interval endpoint here: omit that origin row.
    """
    z = np.asarray(depths, dtype=float)
    V = np.asarray(cumulative_variance, dtype=float)
    if z.ndim != 1 or z.size == 0 or V.shape != z.shape:
        raise ValueError('Depths and variances must be matching nonempty 1D arrays.')
    if (not np.isfinite(R0) or R0 <= 0 or not np.all(np.isfinite(z))
            or np.any(z <= 0) or np.any(np.diff(z) <= 0)):
        raise ValueError('Need positive increasing boundary depths and positive R0.')
    if layer_ids is not None:
        ids = np.asarray(layer_ids)
        if ids.shape != z.shape or not np.array_equal(ids, np.arange(z.size)):
            raise ValueError('Expected one entry per layer with IDs 0,1,...; '
                             'check omitted or duplicated layers.')
    inside = z < R0
    z, V = z[inside].copy(), V[inside].copy()
    V = _variance(V)
    dz = np.diff(np.r_[0.0, z])
    _grid(z, dz, R0)
    return z, dz, V


def highland_variance(z, dz, R0, E0, p_exp, X0):
    """Regrid the ORIGINAL integrated-Highland/global-log prescription.

    This preserves that heuristic prescription; it does not establish its
    validity near stopping or for arbitrary thicknesses. Integrates the
    energy-dependent factor with midpoint quadrature and the actual dz.
    """
    z, dz = _grid(z, dz, R0)
    if not np.all(np.isfinite([E0, p_exp, X0])) or min(E0, p_exp, X0) <= 0:
        raise ValueError('E0, p_exp and X0 must be finite and positive.')
    
    mid = z - dz/2
    energy = _energy_at_depth(mid, R0, E0, p_exp)
        
    if energy.shape != z.shape or np.any(~np.isfinite(energy)) or np.any(energy <= 0):
        raise ValueError('Energy function must return finite positive energies.')
    mass = 938.272
    beta_pc = energy * (energy + 2*mass) / (energy + mass)
    V = (1+0.038*np.log(z/X0))**2 * np.cumsum((13.6/beta_pc)**2 * dz/X0)
    return _variance(V)

def generalized_highland_variance(z, dz, R0, E0, p_exp, X0):
    
    z, dz = _grid(z, dz, R0)
    if not np.all(np.isfinite([E0, p_exp, X0])) or min(E0, p_exp, X0) <= 0:
        raise ValueError('E0, p_exp and X0 must be finite and positive.')
    
    mid = z - dz/2
    energy = _energy_at_depth(mid, R0, E0, p_exp)
    
    if energy.shape != z.shape or np.any(~np.isfinite(energy)) or np.any(energy <= 0):
        raise ValueError('Energy function must return finite positive energies.')
    mass = 938.272
    beta = _beta(mass, energy)
    momentum = _momentum(mass, energy)
    reduced_thickness = dz / (X0*beta**2)
    reduced_totalthickness = z / (X0*beta**2)

    V = (1+0.038*np.log(reduced_totalthickness))**2 * np.cumsum((13.6/momentum)**2 * reduced_thickness)
    return _variance(V)

def scattering_length_from_X0(X0, Z):
    """Fruehwirth-Regler Eq. (37): Xs = X0*h(Z), in the units of X0.

    This is an elemental-material relation. Using an effective Z for a
    compound (e.g. the current PbWO4/water settings) is an approximation.
    """
    if not np.all(np.isfinite([X0, Z])) or min(X0, Z) <= 0:
        raise ValueError('X0 and Z must be finite and positive.')
    h = (Z + 1) / Z * np.log(287 * Z**(-1/2)) / np.log(159 * Z**(-1/3))
    Xs = X0 * h
    if not np.isfinite(Xs) or Xs <= 0:
        raise ValueError('Invalid scattering length from X0 and Z.')
    return float(Xs)


def gaussian_mixture_model(z, dz, R0, E0, p_exp, X0, Z, Xs=None):
    """Fruehwirth-Regler Xs mixture, Eqs. (37), (40), (64), (65).

    Reference: NIM A 456 (2001) 369-389, doi:10.1016/S0168-9002(00)00589-1.
    Returns (per-layer core variance, per-layer tail variance, tail probability).
    Variances are for ONE projected angle in rad^2. Lengths are cm; E0 is MeV.
    The polynomial gives a dimensionless variance, not a standard deviation.

    Xs defaults to X0*h(Z). An explicit Xs may instead be supplied in cm.
    An effective Z for a compound remains an approximation, including in the
    tail-weight fit. The paper's fitted examples are carbon, silicon and iron.
    Energy loss is treated with the existing midpoint-energy approximation
    separately in each layer, not with one constant energy over the full path.
    """
    z, dz = _grid(z, dz, R0)
    if not np.all(np.isfinite([E0, p_exp, X0, Z])) or min(E0, p_exp, X0, Z) <= 0:
        raise ValueError('E0, p_exp, X0 and Z must be finite and positive.')
    Xs = scattering_length_from_X0(X0, Z) if Xs is None else float(Xs)
    if not np.isfinite(Xs) or Xs <= 0:
        raise ValueError('Xs must be finite and positive.')
    
    mid = z - dz/2
    energy = _energy_at_depth(mid, R0, E0, p_exp)
    
    if energy.shape != z.shape or np.any(~np.isfinite(energy)) or np.any(energy <= 0):
        raise ValueError('Energy function must return finite positive energies.')
    mass = 938.272
    beta = _beta(mass, energy)
    
    reduced_thickness = dz / (Xs * beta**2)
    modified_reduced_thickness = Z**(2/3) * reduced_thickness

    log_d = np.log(reduced_thickness)
    log_d_modified = np.log(modified_reduced_thickness)

    variance1 = 0.8510 + 0.03314*log_d - 0.001825*log_d**2

    epsilon = np.where(
        log_d_modified < 0.5,
        0.05283 + 0.0077*log_d_modified + 0.00069*log_d_modified**2,
        -0.01517 + 0.1151*log_d_modified - 0.00653*log_d_modified**2,
    )

    if (
        np.any(~np.isfinite(variance1))
        or np.any(~np.isfinite(epsilon))
        or np.any((variance1 <= 0) | (variance1 >= 1))
        or np.any((epsilon <= 0) | (epsilon >= 0.5))
    ):
        raise ValueError('Xs parametrization gives invalid core/tail parameters; '
                         'check layer thickness, energy and fit applicability.')

    variance2 = (1 - (1 - epsilon)*variance1) / epsilon

    pc = np.sqrt(energy * (energy + 2 * mass))  # MeV

    # Eq. (40): 0.000225 GeV^2 becomes 225 MeV^2 = (15 MeV)^2.
    # No Highland logarithmic correction: this is the full second moment.
    total_variance = (15.0 / (beta * pc))**2 * dz / Xs  # rad²

    variance1 = variance1 * total_variance  # core variance in rad²
    variance2 = variance2 * total_variance  # tail variance in rad²

    return variance1, variance2, epsilon

def weighted_eigenvalues(V, dz):
    """Weights of sum_k w_k chi2_2 for nonuniform integration intervals.

    C_ij = V[min(i,j)] for independent angular increments.
    B = diag(sqrt(dz)) C diag(sqrt(dz)) / 2 is symmetric.
    Its eigenvalues ALREADY include interval widths and the factor 1/2.
    """
    V = _variance(V)
    dz = np.asarray(dz, float)
    if dz.shape != V.shape or np.any(~np.isfinite(dz)) or np.any(dz <= 0):
        raise ValueError('Widths must match V and be finite and positive.')
    indices = np.arange(V.size)
    C = V[np.minimum.outer(indices, indices)]
    root = np.sqrt(dz)
    B = 0.5 * root[:, None] * C * root[None, :]
    weights = np.linalg.eigvalsh(B)
    tolerance = 100*np.finfo(float).eps*V.size*max(np.max(np.abs(weights)),
                                                np.finfo(float).tiny)
    if weights[0] < -tolerance:
        raise ValueError('Weighted covariance has significantly negative eigenvalues.')
    return weights[weights > 0]

def calculate_coefficients(weights):
    """Original residue formula; may be ill-conditioned for nearby weights."""
    weights = np.asarray(weights, dtype=float)
    coefficients = np.ones(len(weights))
    for i, wi in enumerate(weights):
        for j, wj in enumerate(weights):
            if i != j:
                if wi == wj:
                    raise ValueError("Repeated weights: residue CDF is singular.")
                coefficients[i] *= wi / (wi - wj)
    return coefficients

def gchi2_exact_pdf(x, weights):
    weights = np.asarray(weights, dtype=float)
    weights = weights[weights > 0]
    x = np.asarray(x, dtype=float)
    if not weights.size:
        return np.asarray(x >= 0, dtype=float)
    coefficients = calculate_coefficients(weights)
    pdf = np.zeros_like(x)
    for wi, Ai in zip(weights, coefficients):
        # pdf += Ai * np.exp(-np.maximum(x, 0) / (2 * wi)) // non-normalized
        pdf += Ai / (2 * wi) * np.exp(-np.maximum(x, 0) / (2 * wi))
    return np.where(x <= 0, 0.0, pdf)

def gchi2_exact_cdf(x, weights):
    weights = np.asarray(weights, dtype=float)
    weights = weights[weights > 0]
    x = np.asarray(x, dtype=float)
    if not weights.size:
        return np.asarray(x >= 0, dtype=float)
    coefficients = calculate_coefficients(weights)
    cdf = np.ones_like(x)
    
    for wi, Ai in zip(weights, coefficients):
        cdf -= Ai * np.exp(-np.maximum(x, 0) / (2 * wi))
    
    return np.where(x <= 0, 0.0, cdf)

def reach_cdf_grid(z, dz, V, R0, cdf, indices=None):
    """Optional analytical evaluation only at selected grid indices.

    All preceding intervals remain in each covariance matrix even if their
    reach probabilities are not requested. cdf(x, weights) is user-supplied.
    """
    z, dz = _grid(z, dz, R0)
    V = _variance(V)
    if V.shape != z.shape:
        raise ValueError('One cumulative variance is required per grid point.')
    if indices is None:
        indices = np.arange(z.size)
    indices = np.asarray(indices)
    if indices.ndim != 1 or indices.dtype.kind not in 'iu':
        raise ValueError('indices must be a 1D integer array.')
    if np.any(indices < 0) or np.any(indices >= z.size):
        raise ValueError('Evaluation index outside grid.')
    probabilities = np.empty(indices.size)
    for k, j in enumerate(indices):
        w = weighted_eigenvalues(V[:j+1], dz[:j+1])
        probability = 1.0 if w.size == 0 else float(cdf(R0-z[j], w))
        if not np.isfinite(probability) or not 0 <= probability <= 1:
            raise FloatingPointError(f'CDF returned {probability} at {z[j]} cm; '
                                     'check numerical cancellation.')
        probabilities[k] = probability
    return z[indices], probabilities


def _sampling_counts(n_samples, batch_size):
    n_samples, batch_size = operator.index(n_samples), operator.index(batch_size)
    if n_samples <= 0 or batch_size <= 0:
        raise ValueError('Sample and batch counts must be positive integers.')
    return n_samples, batch_size


def gchi2_cdf_rng(x, weights, rng, n_samples=50_000, batch_size=1000):
    """Original per-depth RNG: estimate P(sum_k weights[k]*chi2_2 <= x).

    weights already contain nonuniform interval widths and the factor 1/2.
    Returns (probability, pointwise sampling standard error).
    """
    weights = np.asarray(weights, dtype=float)
    if (weights.ndim != 1 or np.any(~np.isfinite(weights))
            or np.any(weights < 0)):
        raise ValueError('Weights must be a finite nonnegative 1D array.')
    if not np.isfinite(x):
        raise ValueError('Threshold must be finite.')
    n_samples, batch_size = _sampling_counts(n_samples, batch_size)
    weights = weights[weights > 0]
    if x < 0:
        return 0.0, 0.0
    if weights.size == 0:
        return 1.0, 0.0
    if x == 0:
        return 0.0, 0.0
    n_reach = 0
    for start in range(0, n_samples, batch_size):
        size = min(batch_size, n_samples-start)
        samples = rng.chisquare(df=2, size=(size, len(weights)))
        delta_R = samples @ weights
        n_reach += np.count_nonzero(delta_R <= x)
    probability = n_reach/n_samples
    error = np.sqrt(probability*(1-probability)/n_samples)
    return probability, error

def gchi2_cdf_rng_mixture(x, CoreWeights, TailWeights, core_probability, rng, n_samples=50_000, batch_size=1000):
    core = np.asarray(CoreWeights, dtype=float)
    tail = np.asarray(TailWeights, dtype=float)

    for weights in (core, tail):
        if (weights.ndim != 1
                or np.any(~np.isfinite(weights))
                or np.any(weights < 0)):
            raise ValueError("Weights must be finite nonnegative 1D arrays.")

    if not np.isfinite(x):
        raise ValueError("Threshold must be finite.")

    w = float(core_probability)
    if not np.isfinite(w) or not 0 <= w <= 1:
        raise ValueError("core_probability must be between 0 and 1.")

    n_samples, batch_size = _sampling_counts(n_samples, batch_size)
    core = core[core > 0]
    tail = tail[tail > 0]

    if x < 0:
        return 0.0, 0.0

    n_reach = 0

    for start in range(0, n_samples, batch_size):
        size = min(batch_size, n_samples - start)

        # One component choice per complete sample.
        n_core = np.count_nonzero(rng.random(size) < w)
        n_tail = size - n_core

        for count, weights in ((n_core, core), (n_tail, tail)):
            if count == 0:
                continue

            if weights.size == 0:
                n_reach += count
                continue

            samples = rng.chisquare(df=2, size=(count, weights.size))
            delta_R = samples @ weights
            n_reach += np.count_nonzero(delta_R <= x)

    probability = n_reach / n_samples
    error = np.sqrt(probability * (1 - probability) / n_samples)
    return probability, error

def reach_probability_rng_geometry(remaining_range, cumulative_variance, dz,
                                   rng, n_samples=50_000, batch_size=1000):
    """Original per-depth correlated-angle RNG with nonlinear forward geometry.

    Pass a cumulative-variance PREFIX and its corresponding width PREFIX.
    dz may also be a positive scalar for backwards compatibility.
    Returns (probability, pointwise standard error, invalid_fraction).
    A forward-domain violation raises; a successful call returns 0 for the
    invalid_fraction, preserving the original three-value return signature.
    """
    V = _variance(cumulative_variance)
    widths = np.broadcast_to(np.asarray(dz, dtype=float), V.shape)
    if np.any(~np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError('Widths must be finite and positive.')
    if not np.isfinite(remaining_range):
        raise ValueError('Remaining range must be finite.')
    n_samples, batch_size = _sampling_counts(n_samples, batch_size)
    if remaining_range < 0:
        return 0.0, 0.0
    if not np.any(V):
        return 1.0, 0.0
    if remaining_range == 0:
        return 0.0, 0.0
    increment_sigma = np.sqrt(np.diff(np.r_[0.0, V]))
    n_reach = 0
    n_nonforward = 0
    for start in range(0, n_samples, batch_size):
        size = min(batch_size, n_samples-start)
        kicks_x = rng.normal(size=(size, V.size))*increment_sigma
        kicks_y = rng.normal(size=(size, V.size))*increment_sigma
        theta_x = np.cumsum(kicks_x, axis=1)
        theta_y = np.cumsum(kicks_y, axis=1)
        forward = np.all(
            (np.abs(theta_x) < np.pi / 2)
            & (np.abs(theta_y) < np.pi / 2),
            axis=1,
        )

        # Non-forward trajectories count as failures to reach.
        # Evaluate tan() only for trajectories inside its valid domain.
        theta_x = theta_x[forward]
        theta_y = theta_y[forward]

        slope_squared = np.tan(theta_x)**2 + np.tan(theta_y)**2
        excess_factor = slope_squared / (np.sqrt(1 + slope_squared) + 1)
        delta_R = excess_factor @ widths
        
        n_nonforward += np.count_nonzero(~forward)
        n_reach += np.count_nonzero(delta_R <= remaining_range)
    # print(f"Non-forward trajectories: {n_nonforward}")
    probability = n_reach/n_samples
    error = np.sqrt(probability*(1-probability)/n_samples)
    return probability, error

def reach_probability_rng_geometry_mixture(remaining_range, epsilon, cumulative_core_variance, cumulative_tail_variance, dz,
                                   rng, n_samples=50_000, batch_size=1000):
    Vc = _variance(cumulative_core_variance)
    Vt = _variance(cumulative_tail_variance)
    widths = np.broadcast_to(np.asarray(dz, dtype=float), Vc.shape)

    if np.any(~np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError('Widths must be finite and positive.')
    if not np.isfinite(remaining_range):
        raise ValueError('Remaining range must be finite.')
    n_samples, batch_size = _sampling_counts(n_samples, batch_size)
    if remaining_range < 0:
        return 0.0, 0.0
    if not np.any(Vc) and not np.any(Vt):
        return 1.0, 0.0
    if remaining_range == 0:
        return 0.0, 0.0
    
    increment_var_core = np.diff(np.r_[0.0, Vc])
    increment_var_tail = np.diff(np.r_[0.0, Vt])
    n_reach = 0
    n_nonforward = 0
    for start in range(0, n_samples, batch_size):
        size = min(batch_size, n_samples-start)

        shape = (size, len(Vc))

        is_tail = rng.random(shape) < epsilon
        kick_sigma = np.sqrt(np.where(is_tail, increment_var_tail, increment_var_core))

        kicks_x = rng.normal(size=shape) * kick_sigma
        kicks_y = rng.normal(size=shape) * kick_sigma

        theta_x = np.cumsum(kicks_x, axis=1)
        theta_y = np.cumsum(kicks_y, axis=1)
        
        forward = np.all(
            (np.abs(theta_x) < np.pi / 2)
            & (np.abs(theta_y) < np.pi / 2),
            axis=1,
        )

        # Non-forward trajectories count as failures to reach.
        # Evaluate tan() only for trajectories inside its valid domain.
        theta_x = theta_x[forward]
        theta_y = theta_y[forward]

        slope_squared = np.tan(theta_x)**2 + np.tan(theta_y)**2
        excess_factor = slope_squared / (np.sqrt(1 + slope_squared) + 1)
        delta_R = excess_factor @ widths
        
        n_nonforward += np.count_nonzero(~forward)
        n_reach += np.count_nonzero(delta_R <= remaining_range)
    
    probability = n_reach/n_samples
    error = np.sqrt(probability*(1-probability)/n_samples)
    return probability, error

def gchi2_cdf_checked(x, weights):
    weights = np.asarray(weights, dtype=float)
    if (weights.ndim != 1 or np.any(~np.isfinite(weights))
            or np.any(weights < 0) or not np.isfinite(x)):
        raise ValueError('Require finite nonnegative weights and finite x.')
    weights = weights[weights > 0]
    if x < 0:
        return 0.0
    if not weights.size:
        return 1.0
    if x == 0:
        return 0.0
    w = weights.astype(np.longdouble)
    differences = w[:, None] - w[None, :]
    np.fill_diagonal(differences, 1)
    if np.all(differences != 0):
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            ratios = w[:, None] / differences
            np.fill_diagonal(ratios, 1)
            terms = np.prod(ratios, axis=1) * np.exp(-np.longdouble(x)/(2*w))
            probability = 1 - np.sum(terms)
            # Conservative roundoff screen; not a rigorous error bound.
            roundoff = (100 * len(w) * np.finfo(np.longdouble).eps
                        * (1 + np.sum(np.abs(terms))))
        if (np.isfinite(probability) and roundoff < 1e-10
                and 0 <= probability <= 1):
            return float(probability)
    ### secondary approach
    from scipy.linalg import expm
    print("First approach for CDF failed resorting to secondary approach")
    rates = 1 / (2 * weights)
    generator = np.zeros((len(weights) + 1, len(weights) + 1))
    indices = np.arange(len(weights))
    generator[indices, indices] = -rates
    generator[indices, indices + 1] = rates
    probability = float(expm(x * generator)[0, -1])
    if not np.isfinite(probability) or not -1e-10 <= probability <= 1 + 1e-10:
        raise FloatingPointError('Conditional CDF evaluation failed.')
    return float(np.clip(probability, 0, 1))


def reach_probability_sampled_mixture(z, dz, R0, var_core, var_tail,
                                      p_tail, rng, n_configurations=500,
                                      cdf=None, progress=None):
    z, dz = _grid(z, dz, R0)
    var_core, var_tail, p_tail = [np.asarray(a, dtype=float)
                                for a in (var_core, var_tail, p_tail)]
    if any(a.shape != z.shape or np.any(~np.isfinite(a))
           for a in (var_core, var_tail, p_tail)):
        raise ValueError('One finite variance/probability per layer required.')
    if (np.any(var_core < 0) or np.any(var_tail < 0)
            or np.any((p_tail < 0) | (p_tail > 1))):
        raise ValueError('Variances must be nonnegative; probabilities in [0,1].')
    n_configurations = operator.index(n_configurations)
    if n_configurations < 2:
        raise ValueError('Use at least two configurations to estimate an error.')
    cdf = gchi2_cdf_checked if cdf is None else cdf
    choices = rng.random((n_configurations, len(z))) < p_tail
    mean = np.empty(len(z))
    error = np.empty(len(z))
    for j in range(len(z)):
        patterns, counts = np.unique(choices[:, :j+1], axis=0, return_counts=True)
        conditional = np.empty(len(patterns))
        for k, pattern in enumerate(patterns):
            cumulative = np.cumsum(np.where(pattern, var_tail[:j+1], var_core[:j+1]))
            eigenvalues = weighted_eigenvalues(cumulative, dz[:j+1])
            conditional[k] = float(cdf(R0-z[j], eigenvalues))
        if (np.any(~np.isfinite(conditional))
                or np.any((conditional < 0) | (conditional > 1))):
            raise FloatingPointError(f'Invalid conditional CDF at depth {z[j]}.')
        mean[j] = np.dot(counts / n_configurations, conditional)
        error[j] = np.sqrt(np.dot(counts, (conditional - mean[j])**2)
                           / (n_configurations * (n_configurations-1)))
        if progress is not None:
            progress(j+1, len(z), len(patterns))
    return mean, error

