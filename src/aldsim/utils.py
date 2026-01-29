#Copyright © 2024-Present, UChicago Argonne, LLC

import math as m

from .constants import amu, kb, Nav


def calc_vth(M, T):
    """Compute the mean thermal velocity

    Parameters
    ----------

    M : float
        Molecular mass in atomic mass units
    T : float
        Temperature in K

    Returns
    -------
    float
        Mean thermal velocity in m/s

    """
    return m.sqrt(8*kb*T/(m.pi*amu*M))


def sitearea_fromgpc(gpc, M, density, nmol=1):
    """Average area of a surface site from growth per cycle

    Calculate the average area of a surface site

    Parameters
    ----------

    gpc : float
        Growth per cycle, in Angstroms
    M : float
        Molecular mass from the solid in atomic mass units
    density : float
        Density of the film, in g/cm3
    nmol : int, optional (default 1)
        Number of precursor molecules per unit formula of the solid

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    masscm2 = density*gpc*1e-8
    molcm2 = masscm2/M*6.022e23
    return 1e-4/(nmol*molcm2)


def sitearea_fromqcm(mpc, M, nmol=1):
    """Average area of a surface site from QCM mass

    Calculate the average area of a surface site from qcm data

    Parameters
    ----------

    mpc : float
        Mass per cycle in  ng/cm2
    M : float
        Molecular mass in atomic mass units
    nmol : int, optional (default 1)
        Number of precursor molecules per unit formula of the solid

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    return M/(mpc*1e-5*6.022e23*nmol)


def sitedensity_fromqcm(mpc, M, nmol=1):
    """Surface site density from QCM mass

    Calculate the surface site density from qcm data

    Parameters
    ----------

    mpc : float
        Mass per cycle in  ng/cm2
    M : float
        Molar mass in g/mol
    nmol : int, optional (default 1)
        Number of precursor molecules per unit formula of the solid

    Returns
    -------
    float
        Surface site density in #sites / sq. meters

    """

    return 1/sitearea_fromqcm(mpc, M, nmol)


def sitearea_fromrbs(atoms_area, atoms_permol=1.0):
    """Average area of a surface site from RBS data

    Calculate the average area of a surface site from RBS data

    Parameters
    ----------

    atoms_area : float
        Atoms per unit area for one ALD cycle (atoms per sq. meter)
    atoms_permol : int, optional (default 1)
        Number of atoms per precursor molecule

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    return atoms_permol/atoms_area

