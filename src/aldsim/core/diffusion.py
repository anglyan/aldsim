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
    """Model for ALD coating through diffusion in vias or trenches.

    Implementation of a non-dimensional model for atomic layer deposition
    in high-aspect-ratio structures such as vias or trenches. The model
    accounts for precursor diffusion limitations and surface reaction kinetics.

    The model assumes a first-order irreversible Langmuir kinetics
    with the sticking probability value contained in the Damkohler
    number.

    Args:
        AR (float) : Aspect ratio
        p_stick0 (float) : sticking probability
        p_rec (float) : recombination probability

    """

    model_kwd = ["dose", "nondim"]

    def __init__(self, AR, p_stick0, p_rec0=0, p_rec1=0):
        self.AR = AR
        self.p_stick0 = p_stick0
        self.p_rec0 = p_rec0
        self.p_rec1 = p_rec1


    def run(self, tmax=5, dt=0.01):
        """Runs the simulation for a given or predefined amount of time

        Runs the model for a predefined or user-provided time

        Args:
            tmax (float, optional): largest normalized dose time.
            dt (float, optional): time step value.

        Returns:
            A tuple of time, surface coverage, precursor utilization arrays

        """
        pass

    def run_until_cov(self, N=None, max_cov=0.99, save_every=0.2, dt=0.05):
        """Runs the simulation for a given or predefined amount of time

        Runs the model for a predefined or user-provided time

        Args:
            tmax (float, optional): largest normalized dose time.
            dt (float, optional): time step value.

        Returns:
            A tuple of time, surface coverage, precursor utilization arrays

        """        
        if N is None:
            N = int(4*self.AR)
        return solve_until_cov(self.AR, N, self.p_stick0, self.p_rec0, self.p_rec1, max_cov, save_every, dt)



def calc_coverage(Da, t):
    """Analytical expression of the surface coverage for the DiffusionVia model"""
    pass
