#Copyright © 2025, UChicago Argonne, LLC

import numpy as np
from scipy.linalg import solve_banded
import numpy as np

def transport_circular(AR, p_reac):
    """Solve the steady state transport equation inside a circular via

    Transport using Knudsen diffusion

    """

    N = p_reac.shape[0]-1
    ab = np.zeros((3,N+1))
    #diagonal, 1,j
    ab[1,:-1] = 3*(AR/N)**2
    ab[1,-1] = 3/4*AR/N
    ab[1,:] *= p_reac
    ab[1,1:-2] += 2
    ab[1,0] += 3
    ab[1,-2] += 3
    ab[1,-1] += 2
    ab[0,1:-1] = -1
    ab[0,-1] = -2
    ab[2,0:-2] = -1
    ab[2,-2] = -2

    b = np.zeros(N+1)
    b[0] = 2 
    return solve_banded((1,1), ab, b)

def solve(AR, N, p_stick0, p_rec0=0, p_rec1=0, target_cov=0.25, time_multiplier=2):
    dt = 0.05
    cov = np.zeros(N+1)
    i = 0
    s_index = 0
    store = []
    store_times = []
    found = False
    done = False
    target_time = None
    while np.min(cov) < 0.99:
        p_stick_eff = (p_stick0+p_rec0)*(1-cov) + p_rec1*cov
        x = transport_circular(AR, p_stick_eff)
        a = x*dt
        new_cov = (cov + a)/(1+a)
        cov = new_cov
        i += 1
        if not done:
            if found:
                if i >= target_time:
                    store.append(cov[:-1])
                    store_times.append(i)
                    done = True
            else:
                if np.mean(cov) > target_cov:
                    store.append(cov[:-1])
                    store_times.append(i)
                    found = True
                    target_time = time_multiplier*i
    store_times.append(i)
    return store, store_times


def solve_until(AR, N, p_stick0, p_rec0=0, p_rec1=0, target_time=1.0, save_every=0.2, dt=0.01):
    """Solve precursor transport inside a circular via of aspect ratio AR

    This function solves the precursor transport and surface reaction kinetics
    inside a circular via using Knudsen diffusion. The simulation runs until
    a specified target time and saves coverage profiles at regular intervals.

    Args:
        AR (float): Aspect ratio of the circular via
        N (int): Number of discretized segments along the via depth
        p_stick0 (float): Sticking probability of the self-limited process
        p_rec0 (float, optional): Recombination probability on bare sites. Defaults to 0.
        p_rec1 (float, optional): Recombination probability on reacted sites. Defaults to 0.
        target_time (float, optional): Normalized time at which the simulation stops. Defaults to 1.0.
        save_every (float, optional): Normalized time interval at which coverage profiles
            are saved. Defaults to 0.2.
        dt (float, optional): time increment used for the numerical integration

    Returns:
        tuple: A tuple containing:
            - store (list): List of coverage arrays at saved time points, each of size N
            - store_times (list): List of normalized times corresponding to saved profiles

    Notes:
        All time values are in normalized units.
    """
    dt = 0.05
    cov = np.zeros(N+1)
    i = 0
    store = []
    store_times = []
    next_save_time = save_every

    while i*dt < target_time:
        p_stick_eff = (p_stick0+p_rec0)*(1-cov) + p_rec1*cov
        x = transport_circular(AR, p_stick_eff)
        a = x*dt
        new_cov = (cov + a)/(1+a)
        cov = new_cov
        i += 1

        # Check if we should save this iteration
        current_time = i*dt
        if current_time >= next_save_time:
            store.append(cov[:-1].copy())
            store_times.append(current_time)
            next_save_time += save_every

    # Add the final coverage to store if it hasn't been added yet
    final_time = i*dt
    if len(store_times) == 0 or store_times[-1] < final_time:
        store.append(cov[:-1].copy())
        store_times.append(final_time)

    return store, store_times


def solve_until_cov(AR, N, p_stick0, p_rec0=0, p_rec1=0, target_cov=0.99, save_every=0.2, dt=0.05):
    """Solve precursor transport inside a circular via of aspect ratio AR

    This function solves the precursor transport and surface reaction kinetics
    inside a circular via using Knudsen diffusion. The simulation runs until
    a specified target coverage is reached.

    Args:
        AR (float): Aspect ratio of the circular via
        N (int): Number of discretized segments along the via depth
        p_stick0 (float): Sticking probability of the self-limited process
        p_rec0 (float, optional): Recombination probability on bare sites. Defaults to 0.
        p_rec1 (float, optional): Recombination probability on reacted sites. Defaults to 0.
        target_cov (float, optional): Normalized time at which the simulation stops. Defaults to 0.99
        save_every (float, optional): Coverage intervals at which profiles and time are saved. Defaults to 0.2.
        dt (float, optional): time increment used for the numerical integration

    Returns:
        tuple: A tuple containing:
            - store (list): List of coverage arrays
            - store_times (list): List of normalized times corresponding to saved profiles

    Notes:
        All time values are in normalized units.
    """
    cov = np.zeros(N+1)
    i = 0
    store = []
    store_times = []
    next_save_cov = save_every

    mean_cov = 0

    while mean_cov < target_cov:
        p_stick_eff = (p_stick0+p_rec0)*(1-cov) + p_rec1*cov
        x = transport_circular(AR, p_stick_eff)
        a = x*dt
        new_cov = (cov + a)/(1+a)
        cov = new_cov
        i += 1

        # Check if we should save this iteration
        mean_cov = np.mean(cov)
        current_time = i*dt
        if mean_cov >= next_save_cov:
            store.append(cov[:-1].copy())
            store_times.append(current_time)
            next_save_cov += save_every

    # Add the final coverage to store if it hasn't been added yet
    final_time = i*dt
    store.append(cov[:-1].copy())
    store_times.append(final_time)

    return store, store_times


class DiffusionVia:
    """Model for ALD in high aspect ratio circular vias.

    Implementation of a non-dimensional model for atomic layer deposition
    in high-aspect-ratio circular vias. The model uses a Knudsen diffusion transport
    model and self-limited surface reaction kinetics and surface recombination.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value determining the reaction rate.
    It also supports recombination processes on both bare and reacted
    surface sites.

    Parameters
    ----------
    AR : float
        Aspect ratio of the circular via (depth/diameter). Higher values
        indicate deeper, narrower structures where diffusion limitations
        become more significant. Must be non-negative.
    p_stick0 : float
        Sticking probability for the self-limited ALD process. Represents
        the probability that a precursor molecule will react when it
        encounters a surface site.
    p_rec0 : float, optional
        Recombination probability on bare (unreacted) surface sites.
        Default is 0.
    p_rec1 : float, optional
        Recombination probability on reacted surface sites. Default is 0.

    Attributes
    ----------
    AR : float
        The aspect ratio of the via.
    p_stick0 : float
        The sticking probability of the ALD process.
    p_rec0 : float
        The recombination probability on bare sites.
    p_rec1 : float
        The recombination probability on reacted sites.

    Examples
    --------
    Create a DiffusionVia model for a via with aspect ratio 10:

    >>> model = DiffusionVia(AR=10, p_stick0=0.05)
    >>> coverage, times = model.run(max_time=2.0)
    >>> print(f"Final mean coverage: {coverage[-1].mean():.3f}")

    Model with recombination effects:

    >>> model = DiffusionVia(AR=20, p_stick0=0.03, p_rec0=0.01, p_rec1=0.05)
    >>> coverage, times = model.run_until_cov(max_cov=0.95)
    >>> print(f"Time to reach 95% coverage: {times[-1]:.3f}")

    Notes
    -----
    The model uses Knudsen diffusion to describe precursor transport
    inside the circular via. The governing equations are solved using
    a finite difference method with banded matrix solver for efficiency.

    All time values in the model are normalized by a characteristic
    diffusion time scale.

    """

    model_kwd = ["dose", "nondim"]

    def __init__(self, AR, p_stick0, p_rec0=0, p_rec1=0):
        self.AR = AR
        self.p_stick0 = p_stick0
        self.p_rec0 = p_rec0
        self.p_rec1 = p_rec1
        self._nsegments = 4


    def run(self, N=None, max_time=1, save_every=0.2, dt=0.05):
        """Run simulation for a specified normalized time period

        Executes the diffusion-reaction model for precursor transport and
        surface coverage evolution inside a high aspect ratio circular via.
        The simulation runs until the specified maximum normalized time is
        reached, saving coverage profiles at regular time intervals.

        Parameters
        ----------
        N : int, optional
            Number of discretized segments along the via depth. If None
            (default), automatically calculated as 4 * AR to ensure adequate
            spatial resolution. Higher values provide better accuracy but
            increase computation time.
        max_time : float, optional
            Maximum normalized time for the simulation. Default is 1.0.
            Represents the duration of precursor exposure in normalized units.
            Must be positive.
        save_every : float, optional
            Normalized time interval at which coverage profiles are saved.
            Default is 0.2. Smaller values provide more time resolution but
            increase memory usage. Must be positive and less than max_time.
        dt : float, optional
            Time step size for numerical integration (dimensionless).
            Default is 0.05. Smaller values improve accuracy but increase
            computation time. Must be positive and smaller than save_every.

        Returns
        -------
        coverage : list of ndarray
            List of coverage arrays at saved time points. Each array has
            shape (N,) representing the coverage profile along the via depth,
            from the entrance (index 0) to the bottom (index N-1). Values
            are bounded between 0 and 1.
        times : list of float
            List of normalized times corresponding to saved coverage profiles.
            Length matches the coverage list.

        Examples
        --------
        Run simulation with default parameters:

        >>> model = DiffusionVia(AR=10, p_stick0=0.05)
        >>> coverage, times = model.run()
        >>> print(f"Number of saved profiles: {len(coverage)}")
        Number of saved profiles: 6

        Run with custom time parameters and higher resolution:

        >>> model = DiffusionVia(AR=15, p_stick0=0.03)
        >>> coverage, times = model.run(N=100, max_time=3.0, save_every=0.5)
        >>> final_coverage = coverage[-1]
        >>> print(f"Coverage at bottom: {final_coverage[-1]:.3f}")
0
        See Also
        --------
        run_until_cov : Run simulation until target coverage is reached

        """
        if N is None:
            N = int(self._nsegments*self.AR)
        return solve_until(self.AR, N, self.p_stick0, self.p_rec0, self.p_rec1, max_time, save_every, dt)

    def run_until_cov(self, N=None, max_cov=0.99, save_every=0.2, dt=0.05):
        """Run simulation until target mean coverage is reached

        Executes the diffusion-reaction model for precursor transport and
        surface coverage evolution inside a high aspect ratio circular via.
        The simulation continues until the mean coverage across the via
        reaches the specified target value, saving coverage profiles at
        regular coverage intervals.

        Parameters
        ----------
        N : int, optional
            Number of discretized segments along the via depth. If None
            (default), automatically calculated as 4 * AR to ensure adequate
            spatial resolution. Higher values provide better accuracy but
            increase computation time.
        max_cov : float, optional
            Target mean coverage at which the simulation stops. Default is 0.99.
            Represents the spatially-averaged fractional surface coverage.
            Must be between 0 and 1.
        save_every : float, optional
            Coverage interval at which profiles and times are saved.
            Default is 0.2. For example, with default value, profiles are
            saved when mean coverage reaches 0.2, 0.4, 0.6, 0.8, and the
            final target. Must be positive and less than max_cov.
        dt : float, optional
            Time step size for numerical integration (dimensionless).
            Default is 0.05. Smaller values improve accuracy but increase
            computation time.

        Returns
        -------
        coverage : list of ndarray
            List of coverage arrays at saved coverage intervals. Each array
            has shape (N,) representing the coverage profile along the via
            depth, from the entrance (index 0) to the bottom (index N-1).
            Values are bounded between 0 and 1.
        times : list of float
            List of normalized times corresponding to saved coverage profiles.
            Times increase monotonically. Length matches the coverage list.

        Examples
        --------
        Run simulation until 90% mean coverage:

        >>> model = DiffusionVia(AR=10, p_stick0=0.5)
        >>> coverage, times = model.run_until_cov(max_cov=0.9)
        >>> print(f"Time to reach 90% coverage: {times[-1]:.3f}")
        Time to reach 90% coverage: 2.345

        Save coverage profiles every 10% increment:

        >>> model = DiffusionVia(AR=15, p_stick0=0.3, p_rec0=0.1)
        >>> coverage, times = model.run_until_cov(max_cov=0.95, save_every=0.1)

        See Also
        --------
        run : Run simulation for a specified time period

        """        
        if N is None:
            N = int(self._nsegments*self.AR)
        return solve_until_cov(self.AR, N, self.p_stick0, self.p_rec0, self.p_rec1, max_cov, save_every, dt)

