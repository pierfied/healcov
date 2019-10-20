import numpy as np
import healpy as hp
import subprocess
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


def build_cov(cl, nside, mask=None, tree_depth=0, lmax=None, apply_pixwin=False, ninterp=10000, log=False, shift=None):
    tree_depth = 2 ** (tree_depth)

    if mask is None:
        mask = np.ones(hp.nside2npix(nside), dtype=bool)

    if lmax is not None and lmax > len(cl) - 1:
        lmax = len(cl) - 1

        cl = cl[:lmax + 1]

    if apply_pixwin:
        pw = hp.pixwin(nside * tree_depth, lmax=lmax)

        cl = cl[:len(pw)] * (pw ** 2)

    npix = int(mask.sum())
    mask_inds = hp.ring2nest(nside, np.arange(hp.nside2npix(nside))[mask])

    thetas = np.linspace(0, np.pi, ninterp)
    xis = cl2xi_theta(cl, thetas)

    exe_path = os.path.join(os.path.dirname(__file__), 'healcov')
    proc = subprocess.Popen(exe_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    proc.stdin.write(nside.to_bytes(8, 'little'))
    proc.stdin.write(npix.to_bytes(8, 'little'))
    proc.stdin.write(mask_inds.tobytes())
    proc.stdin.write(tree_depth.to_bytes(8, 'little'))
    proc.stdin.write(ninterp.to_bytes(8, 'little'))
    proc.stdin.write(thetas.tobytes())
    proc.stdin.write(xis.tobytes())
    proc.stdin.close()

    cov = np.frombuffer(proc.stdout.read()).reshape([npix, npix])

    info = proc.wait()

    if info != 0:
        raise Exception()

    if log is True:
        cov = np.log(cov / (shift ** 2) + 1)

    return cov
