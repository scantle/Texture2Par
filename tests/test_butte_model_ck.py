"""
test_butte_model_ck.py
=====================
Tests for Texture2Par's cokriging behaviour using the Butte model dataset.

Section 1 – OK reference tests (mirrors test_butte_model.py; kept here so
            that CK-mode regressions can be compared against a known-good
            ordinary-kriging baseline in the same file).
Section 2 – CK-specific tests using the secondary AEM resistivity dataset.
"""

from io import StringIO
import math
import pytest
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.spatial.transform import Rotation

# ---------------------------------------------------------------------------
# Variogram defaults
# ---------------------------------------------------------------------------

_VGM_sph = "sph 0.0 0.12 5000 5000 0.0 100"

# LMC variograms for cokriging (single nested spherical structure).
# Valid LMC: |b_12| = 0.048 ≤ sqrt(b_11 × b_22) = sqrt(0.12 × 0.12) = 0.12
_VGM_P = "sph 0.0 0.12  5000 5000 0.0 150"   # primary  – Coarse fraction
_VGM_S = "sph 0.0 0.12  5000 5000 0.0 150"   # secondary – log10(resistivity)
_VGM_X = "sph 0.0 0.048 5000 5000 0.0"        # cross     – no nnear token required

_covfunc = dict(
    lin = lambda hr: np.where(hr < 1.0, 1.0 - hr, 0.0),
    exp = lambda hr: np.exp(-3.0 * hr),
    log = lambda hr: np.where(hr < 1.0, 1.0 - np.log(1 + hr*(math.e-1.0)), 0.0),
    gau = lambda hr: np.exp(-3.0625 * hr*hr),
    qua = lambda hr: np.where(hr < 1.0, 1.0 - (2.0-hr)*hr, 0.0),
    sph = lambda hr: np.where(hr < 1.0, 1.0 - hr*(1.5-.5*hr*hr), 0.0),
    hol = lambda hr: (1.0 - 3.0 * hr) * np.exp(-3.0 * hr),
    pow = lambda hr: np.where(hr < 1.0, 1.0 - hr**1.5, 0.0),
)


def parse_vgm(svgm, is3d=False):
    """
    Parse a variogram string into components.

    Returns
    -------
    vtype  : str   – covariance model key (3-char lower-case)
    nugget : float
    sill   : float
    rot    : Rotation – anisotropy rotation object
    hh     : ndarray (1, 3) – [hmin, hmax, hvert] range scaling vector
    nmax   : int   – nearest-neighbour limit (9999 if not present in string)
    """
    v = svgm.split()
    vtype = v[0].lower()[:3]
    if vtype not in _covfunc:
        raise ValueError(f"Unknown covariance function {vtype!r}")
    im = 0
    nugget  = float(v[im:=im+1])
    sill    = float(v[im:=im+1])
    hmin    = float(v[im:=im+1])
    hmax    = float(v[im:=im+1])
    if is3d:
        hvert = float(v[im:=im+1])
    else:
        hvert = hmin
    azimuth = float(v[im:=im+1])
    if is3d:
        dip    = float(v[im:=im+1])
        plunge = float(v[im:=im+1])
    else:
        dip = plunge = 0.0
    try:
        nmax = int(v[im:=im+1])
    except (IndexError, ValueError):
        nmax = 9999                   # cross-variograms have no nnear token
    rot = Rotation.from_euler('zyx', [azimuth, dip, plunge], degrees=True)
    hh  = np.array([[hmin, hmax, hvert]])
    return vtype, nugget, sill, rot, hh, nmax


# ---------------------------------------------------------------------------
# Python reference: ordinary kriging (Cholesky Schur complement)
# ---------------------------------------------------------------------------

def ok(obsloc, obsval, newloc, vgm=_VGM_sph, is3d=False):
    """
    OK reference using Cholesky Schur complement.

    Matches T2P for well-conditioned systems; diverges on near-singular
    matrices (e.g. Gaussian variogram with no nugget and dense wells) because
    T2P uses an unpivoted LDL^T decomposition.
    """
    obsloc = np.asarray(obsloc)
    vtype, nugget, sill, rot, hh, nmax = parse_vgm(vgm, is3d=is3d)
    cov0 = nugget + sill

    if not is3d:
        obsloc = np.column_stack([obsloc, np.zeros(len(obsloc))])
        newloc = np.column_stack([newloc, np.zeros(len(newloc))])
    obsloc = rot.apply(obsloc) / hh
    newloc = rot.apply(newloc) / hh

    nobs = len(obsloc)
    if nmax < nobs:
        raise ValueError(f"nmax ({nmax}) must be >= nobs ({nobs})")

    C = _covfunc[vtype](squareform(pdist(obsloc))) * sill
    np.fill_diagonal(C, cov0)
    b = _covfunc[vtype](cdist(obsloc, newloc)) * sill

    cf     = cho_factor(C, lower=True)
    Cinv_b = cho_solve(cf, b)
    Cinv_e = cho_solve(cf, np.ones(nobs))
    mu     = (np.ones(nobs) @ Cinv_b - 1.0) / (np.ones(nobs) @ Cinv_e)
    w      = Cinv_b - np.outer(Cinv_e, mu)

    est = obsval @ w
    var = np.maximum(0.0, cov0 - (w * b).sum(axis=0) - mu)
    return est, var


# ---------------------------------------------------------------------------
# Python reference: ordinary cokriging (single unbiasedness constraint)
# ---------------------------------------------------------------------------

def ck(obsloc0, obsloc1, obsval0, obsval1,
       newloc, vgm0=_VGM_P, vgm1=_VGM_S, vgmc=_VGM_X, is3d=False):
    """
    2-D ordinary cokriging with a single unbiasedness constraint
    (Σ primary weights + Σ secondary weights = 1), matching T2P's
    implementation in kriging.f90 (pset_cok / predict_point_ok).

    The augmented CK system solved for each target x₀ is::

        [C11  C12  1 ] [w1]   [c10(x₀)]
        [C21  C22  1 ] [w2] = [c20(x₀)]
        [1ᵀ   1ᵀ   0 ] [μ ]   [1      ]

    where C11/C22 are direct covariances (vgm0/vgm1), C12 is the cross-
    covariance (vgmc), c10 uses vgm0, c20 uses vgmc, and μ is the Lagrange
    multiplier.  Variance: σ²(x₀) = C11(0) − wᵀ·b − μ.

    All three variograms must share the same rotation and range parameters
    (linear model of coregionalisation with one nested structure).

    Parameters
    ----------
    obsloc0, obsval0 : primary observations  – (n1, 2) and (n1,)
    obsloc1, obsval1 : secondary observations – (n2, 2) and (n2,)
    newloc           : target locations        – (nnew, 2)
    vgm0             : direct variogram string for primary variable
    vgm1             : direct variogram string for secondary variable
    vgmc             : cross-variogram string  (nnear token optional)

    Returns
    -------
    est : (nnew,) cokriging estimates
    var : (nnew,) cokriging variances (clipped to ≥ 0)
    """
    obsloc0 = np.asarray(obsloc0)
    obsloc1 = np.asarray(obsloc1)

    vtype0, nugget0, sill0, rot0, hh0, nmax0 = parse_vgm(vgm0, is3d=is3d)
    vtype1, nugget1, sill1, rot1, hh1, nmax1 = parse_vgm(vgm1, is3d=is3d)
    vtypec, nuggetc, sillc, rotc, hhc, nmaxc = parse_vgm(vgmc, is3d=is3d)

    # All three variograms must share the same spatial structure
    assert np.allclose(rot0.as_matrix(), rot1.as_matrix(), atol=1e-6), \
        "All variograms must share the same rotation"
    assert np.allclose(rot0.as_matrix(), rotc.as_matrix(), atol=1e-6), \
        "All variograms must share the same rotation"
    assert np.allclose(hh0, hh1, atol=1e-6), \
        "All variograms must share the same range parameters"
    assert np.allclose(hh0, hhc, atol=1e-6), \
        "All variograms must share the same range parameters"

    if not is3d:
        obsloc0 = np.column_stack([obsloc0, np.zeros(len(obsloc0))])
        obsloc1 = np.column_stack([obsloc1, np.zeros(len(obsloc1))])
        newloc  = np.column_stack([np.asarray(newloc), np.zeros(len(newloc))])

    # Apply shared rotation and range scaling
    obsloc0 = rot0.apply(obsloc0) / hh0
    obsloc1 = rot0.apply(obsloc1) / hh0
    newloc  = rot0.apply(newloc)  / hh0

    n1, n2, nnew = len(obsloc0), len(obsloc1), len(newloc)
    ntot = n1 + n2 + 1          # +1 for single Lagrange multiplier

    # ---- obs-obs covariance blocks ----------------------------------------
    cov0_p = nugget0 + sill0    # primary sill (used in variance formula)
    C11 = _covfunc[vtype0](squareform(pdist(obsloc0))) * sill0
    np.fill_diagonal(C11, cov0_p)

    C22 = _covfunc[vtype1](squareform(pdist(obsloc1))) * sill1
    np.fill_diagonal(C22, nugget1 + sill1)

    C12 = _covfunc[vtypec](cdist(obsloc0, obsloc1)) * sillc  # (n1, n2), no diag fill

    # ---- augmented CK matrix: [C11  C12  1] / [C21  C22  1] / [1ᵀ  1ᵀ  0] ---
    A = np.zeros((ntot, ntot))
    A[:n1,       :n1]       = C11
    A[:n1,       n1:n1+n2]  = C12
    A[n1:n1+n2,  :n1]       = C12.T
    A[n1:n1+n2,  n1:n1+n2]  = C22
    A[:n1+n2,    ntot-1]    = 1.0
    A[ntot-1,    :n1+n2]    = 1.0

    # ---- RHS: primary-target uses direct vgm; secondary-target uses cross vgm -
    c10 = _covfunc[vtype0](cdist(obsloc0, newloc)) * sill0   # (n1, nnew)
    c20 = _covfunc[vtypec](cdist(obsloc1, newloc)) * sillc   # (n2, nnew)

    B = np.zeros((ntot, nnew))
    B[:n1]      = c10
    B[n1:n1+n2] = c20
    B[ntot-1]   = 1.0

    # ---- solve CK system ---------------------------------------------------
    w = np.linalg.solve(A, B)      # (ntot, nnew)


    # ---- estimates and variance --------------------------------------------
    # shifts secondary values so they share the primary mean before
    # multiplying by weights (Interpolators.f90, line ~221):
    est = obsval0 @ w[:n1] + (obsval1 + obsval0.mean() - obsval1.mean()) @ w[n1:n1+n2]
    var = np.maximum(0.0,
                     cov0_p
                     - (w[:n1+n2] * B[:n1+n2]).sum(axis=0)
                     - w[ntot-1])
    return est, var


# ---------------------------------------------------------------------------
# T2P input-file helpers
# ---------------------------------------------------------------------------

def _t2p_ok(options_extra="", classes="Coarse\n  Fine", datafile="pc2d.csv",
             gridfile="grid.grd", variogram=_VGM_sph) -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_VARIANCE
          INFER_LAST_CLASS
          {options_extra}
        END OPTIONS

        BEGIN CLASSES
          {classes}
        END CLASSES

        BEGIN DATASET
          FILE {datafile}
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE {gridfile}
          XOFFSET   584646
          YOFFSET   4384029.7
          ROTATION  35
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS {classes.split()[0]}
            1 {variogram}
        END VARIOGRAMS
        """


def _t2p_ck(classes="Coarse\n  Fine", covariate="Rho",
             datafile_p="pc2d.csv", datafile_s="aem2d.csv",
             gridfile="grid.grd",
             vgm_p=_VGM_P, vgm_s=_VGM_S, vgm_x=_VGM_X,
             transform_s="LOG10",
             max_sdist=None,
             options_extra="") -> str:
    primary_class = classes.split()[0]
    sdist_line = f"MAX_SECONDARY_DIST {max_sdist}" if max_sdist is not None else ""
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_VARIANCE
          INFER_LAST_CLASS
          {sdist_line}
          {options_extra}
        END OPTIONS

        BEGIN CLASSES
          {classes}
        END CLASSES

        BEGIN COVARIATES
          {primary_class} {covariate}
        END COVARIATES

        BEGIN DATASET
          FILE {datafile_p}
        END DATASET

        BEGIN DATASET
          FILE {datafile_s}
          TRANSFORM {transform_s}
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE {gridfile}
          XOFFSET   584646
          YOFFSET   4384029.7
          ROTATION  35
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS {primary_class}
            1 {vgm_p}
          CLASS {covariate}
            1 {vgm_s}
          CLASS {primary_class}:{covariate}
            1 {vgm_x}
        END VARIOGRAMS
        """


def calc_ck_correlation(run_t2p, butte_ck_data, vgm_p=_VGM_P, vgm_s=_VGM_S, vgm_x=_VGM_X,):
    """
    CK estimates and variances correlate > 0.99 with a Python reference.

    A reproducible 150-point subsample of the AEM data is used for both
    T2P and the Python reference.  With nnear_s = 200 > 150 and
    MAX_SECONDARY_DIST = 1e10, T2P uses every secondary point in the
    subsample for every target — exactly matching the Python reference.

    The CK system solved here is 195 × 195 (44 primary + 150 secondary
    + 1 Lagrange), which is both tractable and well-conditioned for a
    spherical variogram.
    """
    # Reproducible 150-point subsample of the AEM data
    aem_all = pd.read_csv(StringIO(butte_ck_data["files"]["aem2d.csv"]))
    rng     = np.random.default_rng(42)
    idx     = rng.choice(len(aem_all), 150, replace=False)
    aem_sub = aem_all.iloc[idx].reset_index(drop=True)
    aem_sub["ID"] = np.arange(1, len(aem_sub) + 1)
    extra = {**butte_ck_data["files"], "aem2d.csv": aem_sub.to_csv(index=False, lineterminator='\n')}
    r = run_t2p(
        _t2p_ck(vgm_s=vgm_s, vgm_x=vgm_x, vgm_p=vgm_p, max_sdist=1e10),
        extra_files=extra,
    )

    obs = pd.read_csv(StringIO(butte_ck_data["files"]["pc2d.csv"]))
    t2p = pd.read_csv(r.workdir / "test_Coarse.csv")

    # T2P shifts secondary values so they share the primary mean before
    # multiplying by weights (Interpolators.f90, line ~221):
    #   v_adjust = mean(primary) - mean(secondary)
    # Apply the same adjustment here so the Python estimate matches T2P.
    s_log10 = np.log10(aem_sub["Rho"].values)
    v_adjust = obs["Coarse"].mean() - s_log10.mean()
    s_val = s_log10 + v_adjust
    est, var = ck(
        obs[["x", "y"]].values, aem_sub[["x", "y"]].values,
        obs["Coarse"].values,   s_val,
        t2p[["X", "Y"]].values,
    )
    est = np.clip(est, 0, 1)
    corr_est = np.corrcoef(est, t2p["Layer1"].values)[0, 1]
    corr_var = np.corrcoef(var, t2p["Var_Layer1"].values)[0, 1]
    return corr_est, corr_var

# ---------------------------------------------------------------------------
# 1. Cokriging tests (primary + secondary AEM dataset)
# ---------------------------------------------------------------------------

class TestCokrigingButteCK:

    def test_ck_produces_output(self, run_t2p, butte_ck_data):
        """CK runs without error and the expected output CSV is produced."""
        r = run_t2p(_t2p_ck(), extra_files=butte_ck_data["files"])
        assert r.returncode == 0, \
            f"T2P exited with code {r.returncode}:\n{r.combined}"
        assert r.output_exists("test_Coarse.csv"), \
            "Output file test_Coarse.csv not found"

    def test_ck_estimates_in_range(self, run_t2p, butte_ck_data):
        """All CK Coarse fraction estimates lie in [0, 1]."""
        r = run_t2p(_t2p_ck(), extra_files=butte_ck_data["files"])
        t2p = pd.read_csv(r.workdir / "test_Coarse.csv")
        lo, hi = t2p["Layer1"].min(), t2p["Layer1"].max()
        assert lo >= 0.0 and hi <= 1.0, \
            f"CK estimates outside [0, 1]: min={lo:.4f}, max={hi:.4f}"

    def test_ck_variance_nonnegative(self, run_t2p, butte_ck_data):
        """CK estimation variance is non-negative at every grid node."""
        r = run_t2p(_t2p_ck(), extra_files=butte_ck_data["files"])
        t2p = pd.read_csv(r.workdir / "test_Coarse.csv")
        neg = (t2p["Var_Layer1"] < 0).sum()
        assert neg == 0, f"{neg} grid nodes have negative CK variance"

    def test_ck_differs_from_ok(self, run_t2p, butte_ok_data, butte_ck_data):
        """
        CK estimates differ meaningfully from OK estimates.

        AEM resistivity correlates with Coarse fraction (r ≈ 0.40 at nearest
        AEM sample), so the secondary data shifts the kriging solution.
        The mean absolute deviation between CK and OK must exceed 0.005.
        """
        r_ok = run_t2p(_t2p_ok(**butte_ok_data["kws"]), extra_files=butte_ok_data["files"])
        ok_est = pd.read_csv(r_ok.workdir / "test_Coarse.csv")["Layer1"].values
        r_ck = run_t2p(_t2p_ck(), extra_files=butte_ck_data["files"])
        ck_est = pd.read_csv(r_ck.workdir / "test_Coarse.csv")["Layer1"].values
        mad = np.mean(np.abs(ck_est - ok_est))
        assert mad > 0.005, (
            f"CK and OK estimates are nearly identical (MAD={mad:.5f}); "
            "secondary data appears to have no effect"
        )

    def test_ck_variance_less_than_ok(self, run_t2p, butte_ok_data, butte_ck_data):
        """
        Mean CK variance is lower than mean OK variance.

        The dense AEM dataset provides secondary information at virtually
        every grid node; cokriging cannot increase average uncertainty.
        """
        r_ok = run_t2p(_t2p_ok(**butte_ok_data["kws"]), extra_files=butte_ok_data["files"])
        ok_var = pd.read_csv(r_ok.workdir / "test_Coarse.csv")["Var_Layer1"].values
        r_ck = run_t2p(_t2p_ck(), extra_files=butte_ck_data["files"])
        ck_var = pd.read_csv(r_ck.workdir / "test_Coarse.csv")["Var_Layer1"].values
        assert ck_var.mean() < ok_var.mean(), (
            f"CK mean variance {ck_var.mean():.5f} ≥ OK mean variance "
            f"{ok_var.mean():.5f}"
        )

    def test_ck_reference_correlation(self, run_t2p, butte_ck_data):
        """
        CK estimates and variances correlate > 0.99 with a Python reference.
        """

        corr_est, corr_var = calc_ck_correlation(run_t2p, butte_ck_data)
        assert corr_est > 0.99, (
            f"CK estimate correlation with reference = {corr_est:.4f} (expected > 0.99)"
        )
        assert corr_var > 0.99, (
            f"CK variance correlation with reference = {corr_var:.4f} (expected > 0.99)"
        )

    def test_ck_reference_correlation_aniso(self, run_t2p, butte_ck_data):
        """
        CK estimates and variances correlate > 0.99 with a Python reference with anisotropic variograms.
        """
        vgm_p=_VGM_P.replace("5000 5000 0.0", "1000 5000 40.0")
        vgm_s=_VGM_S.replace("5000 5000 0.0", "1000 5000 40.0")
        vgm_x=_VGM_X.replace("5000 5000 0.0", "1000 5000 40.0")
        corr_est, corr_var = calc_ck_correlation(run_t2p, butte_ck_data, vgm_p, vgm_s, vgm_x)
        assert corr_est > 0.99, (
            f"CK estimate correlation with reference = {corr_est:.4f} (expected > 0.99)"
        )
        assert corr_var > 0.99, (
            f"CK variance correlation with reference = {corr_var:.4f} (expected > 0.99)"
        )
