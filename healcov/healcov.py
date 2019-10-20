import numpy as np
import healpy as hp
import ctypes
import os


def cl2xi_theta(cl, theta):
    """Angular correlation function at separation theta from power spectrum.
    Computes the covariance between pixels at separation of theta from provided power spectrum.
    See https://arxiv.org/pdf/1602.08503.pdf equation 20.
    :param cl: Power spectrum.
    :type cl: array-like (float)
    :param theta: Separation angle in radians.
    :type theta: float or array-like (float)
    :return: xi(theta) - Angular correlation at separation theta.
    :rtype: array-like (float)
    """

    # Convert all array-like input to ndarrays.
    cl = np.asarray(cl)
    theta = np.asarray(theta)

    # Check input sizes.
    if cl.ndim != 1:
        raise Exception('Cl must be a 1D array.')
    if theta.ndim > 1:
        raise Exception('Theta must be a 0D or 1D array.')

    # Get array of l values.
    ells = np.arange(0, len(cl))

    # Compute xi(theta) using Legendre polynomials.
    xi = 1 / (4 * np.pi) * np.polynomial.legendre.legval(np.cos(theta), (2 * ells + 1) * cl)

    return xi


class HealcovArgs(ctypes.Structure):
    _fields_ = [
        ('nside', ctypes.c_long),
        ('npix', ctypes.c_long),
        ('mask_inds', ctypes.POINTER(ctypes.c_long)),
        ('tree_depth', ctypes.c_long),
        ('nsamps', ctypes.c_long),
        ('theta_samps', ctypes.POINTER(ctypes.c_double)),
        ('xi_samps', ctypes.POINTER(ctypes.c_double)),
    ]


def build_cov(cl, nside, mask, tree_depth=0, lmax=None, apply_pixwin=False, ninterp=10000, log=False, shift=None):
    tree_depth = 2 ** (tree_depth)

    if lmax is not None and lmax > len(cl) - 1:
        lmax = len(cl) - 1

        cl = cl[:lmax + 1]

    if apply_pixwin:
        pw = hp.pixwin(nside * tree_depth, lmax=lmax)

        cl = cl[:len(pw)] * (pw ** 2)

    npix = mask.sum()
    mask_inds = np.arange(hp.nside2npix(nside))[mask]

    thetas = np.linspace(0, np.pi, ninterp)
    xis = cl2xi_theta(cl, thetas)

    lptr = ctypes.POINTER(ctypes.c_long)
    dptr = ctypes.POINTER(ctypes.c_double)

    args = HealcovArgs()
    args.nside = nside
    args.npix = npix
    args.mask_inds = mask_inds.ctypes.data_as(lptr)
    args.tree_depth = tree_depth
    args.nsamps = ninterp
    args.theta_samps = thetas.ctypes.data_as(dptr)
    args.xi_samps = xis.ctypes.data_as(dptr)

    lib_path = os.path.join(os.path.dirname(__file__), 'libhealcov.so')
    lib = ctypes.cdll.LoadLibrary(lib_path)

    cpp_build_cov = lib.build_cov()
    cpp_build_cov.argtypes = [HealcovArgs]
    cpp_build_cov.restype = dptr

    cpp_res = cpp_build_cov(args)

    cov = np.ctypeslib.as_array(cpp_res, shape=(npix, npix))

    if log is True:
        cov = np.log(cov / (shift ** 2) + 1)

    return cov
