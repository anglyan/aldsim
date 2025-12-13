from aldsim.core.diffusion import transport_circular, solve, solve_until
import numpy as np
import pytest


def test_transport_circular():
    AR = 10
    p_reac = np.full(20, 0.01)
    result = transport_circular(AR, p_reac)
    # Test that function executes without errors
    assert result is not None
    # Test that result is a numpy array
    assert isinstance(result, np.ndarray)
    # Test that result has the correct shape
    assert result.shape == p_reac.shape


def test_solve():
    AR = 20
    N = 40
    p_stick0 = 1e-3
    store, store_times = solve(AR, N, p_stick0)
    # Test that function executes without errors
    assert store is not None
    assert store_times is not None
    # Test that results are lists
    assert isinstance(store, list)
    assert isinstance(store_times, list)
    # Test that store_times is not empty and contains positive values
    assert len(store_times) > 0
    assert all(t > 0 for t in store_times)
    # Test that store_times is monotonically increasing
    assert all(store_times[i] < store_times[i+1] for i in range(len(store_times)-1))
    # Test that store contains arrays
    if len(store) > 0:
        assert all(isinstance(arr, np.ndarray) for arr in store)
        # Test that stored arrays have correct size (N, not N+1)
        assert all(arr.shape[0] == N for arr in store)


def test_solve_until():
    AR = 20
    N = 40
    p_stick0 = 1e-3
    target_time = 1.0
    save_every = 0.2
    store, store_times = solve_until(AR, N, p_stick0, target_time=target_time, save_every=save_every)
    # Test that function executes without errors
    assert store is not None
    assert store_times is not None
    # Test that results are lists
    assert isinstance(store, list)
    assert isinstance(store_times, list)
    # Test that store_times is not empty
    assert len(store_times) > 0
    # Test that all times are less than or equal to target_time
    assert all(t <= target_time for t in store_times)
    # Test that store_times is monotonically increasing
    assert all(store_times[i] < store_times[i+1] for i in range(len(store_times)-1))
    # Test that store contains arrays
    assert all(isinstance(arr, np.ndarray) for arr in store)
    # Test that stored arrays have correct size (N, not N+1)
    assert all(arr.shape[0] == N for arr in store)
    # Test that saves happen approximately at multiples of save_every
    # (allowing for small differences due to discrete time steps)
    dt = 0.05  # same as in the function
    for i in range(len(store_times)-1):  # exclude the last one which is the final save
        # Check that time is close to a multiple of save_every
        expected_multiple = round(store_times[i] / save_every)
        assert abs(store_times[i] - expected_multiple * save_every) < dt
