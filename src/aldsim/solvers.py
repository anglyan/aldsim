#Copyright © 2024-Present, UChicago Argonne, LLC

from scipy.integrate import solve_ivp

def ode_solver(fdot, initial, tmax, t_eval):
    """
    Wrapper for ODE solver from scipy.integrate
    """
    return solve_ivp(fdot, [0,tmax], initial, t_eval=t_eval, method='LSODA')


def boundedNewton_solver(f, fdot):
    """
    Solves a nonlinear equation bounded between 0 and 1
    """
    ep = 1
    t = 0.5
    damp = 0.1
    while ep > 1e-6:
        f_t = f(t)
        fp_t = fdot(t)
        damp = 0.1
        tn = 1.1
        while tn > 1 or tn < 0:
            damp = 0.5*damp
            tn = t - damp*f_t/fp_t
        ep = abs(t-tn)/t
        t = tn
    return t