"""
test_05_butte_model_ck.py
==========================
Tests for Texture2Par's cokriging behaviour using the Butte model dataset.

Section 1 – OK reference tests (mirrors test_04_butte_model.py; kept here so
            that CK-mode regressions can be compared against a known-good
            ordinary-kriging baseline in the same file).
Section 2 – CK-specific tests using the secondary AEM resistivity dataset.
Section 3 – CK-specific tests extended to 3D space with random Z coordinates.
"""

from io import StringIO
import math
import pytest
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.spatial.transform import Rotation
import time

# ---------------------------------------------------------------------------
# Variogram defaults
# ---------------------------------------------------------------------------

_VGM_sph = "sph 0.0 0.12 5000 5000 0.0 100"

# LMC variograms for cokriging (single nested spherical structure).
# Valid LMC: |b_12| = 0.048 ≤ sqrt(b_11 × b_22) = sqrt(0.12 × 0.12) = 0.12
_VGM_P = "sph 0.0 0.12  5000 5000 0.0 150"   # primary  – Coarse fraction
_VGM_S = "sph 0.0 0.12  5000 5000 0.0 150"   # secondary – log10(resistivity)
_VGM_X = "sph 0.0 0.048 5000 5000 0.0"        # cross     – no nnear token required

# ---------------------------------------------------------------------------
# Variogram defaults (3D Extension)
# ---------------------------------------------------------------------------

# LMC variograms for 3D cokriging (Format: vtype nugget sill hmin hmax hvert azimuth dip plunge [nmax])
# Valid LMC: |b_12| = 0.048 <= sqrt(b_11 * b_22) = 0.12
_VGM_P_3D = "sph 0.0 0.12  5000 5000 1500 35.0 10.0 5.0 150"   # primary  – Coarse fraction
_VGM_S_3D = "sph 0.0 0.12  5000 5000 1500 35.0 10.0 5.0 150"   # secondary – log10(resistivity)
_VGM_X_3D = "sph 0.0 0.048 5000 5000 1500 35.0 10.0 5.0"        # cross

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
    nugget  = float(v[(im:=im+1)])
    sill    = float(v[(im:=im+1)])
    hmin    = float(v[(im:=im+1)])
    hmax    = float(v[(im:=im+1)])
    if is3d:
        hvert = float(v[(im:=im+1)])
    else:
        hvert = hmin
    azimuth = float(v[(im:=im+1)])
    if is3d:
        dip    = float(v[(im:=im+1)])
        plunge = float(v[(im:=im+1)])
    else:
        dip = plunge = 0.0
    try:
        nmax = int(v[(im:=im+1)])
    except (IndexError, ValueError):
        nmax = 9999                   # cross-variograms have no nnear token
    rot = Rotation.from_euler('zxy', [azimuth, -dip, -plunge], degrees=True)
    hh  = np.array([[hmin, hmax, hvert]])
    return vtype, nugget, sill, rot, hh, nmax


def _parse_vgm_ranges_angles(svgm, is3d=False):
    v = svgm.split()
    hmin = float(v[3])
    hmax = float(v[4])
    if is3d:
        hvert = float(v[5])
        azimuth = float(v[6])
        dip = float(v[7])
        plunge = float(v[8])
    else:
        hvert = hmin
        azimuth = float(v[5])
        dip = plunge = 0.0
    return hmin, hmax, hvert, azimuth, dip, plunge


def _t2p_variogram_rotmat(hmin, hmax, hvert, azimuth, dip, plunge):
    # Mirrors T2P's current Fortran anisotropy convention.
    rotmat = np.diag([hmax / hmin, 1.0, hmax / hvert])

    if abs(plunge) > 1e-10:
        alpha = np.deg2rad(-plunge)
        sina = np.sin(alpha)
        cosa = np.cos(alpha)
        tmp = np.array([
            [cosa, 0.0, sina],
            [0.0, 1.0, 0.0],
            [-sina, 0.0, cosa],
        ])
        rotmat = rotmat @ tmp

    if abs(dip) > 1e-10:
        alpha = np.deg2rad(dip)
        sina = np.sin(alpha)
        cosa = np.cos(alpha)
        # Positive dip follows T2P's current positive-downward Fortran
        # matrix convention.
        tmp = np.array([
            [1.0, 0.0, 0.0],
            [0.0, cosa, -sina],
            [0.0, sina, cosa],
        ])
        rotmat = rotmat @ tmp

    if abs(azimuth) > 1e-10:
        alpha = np.deg2rad(azimuth)
        sina = np.sin(alpha)
        cosa = np.cos(alpha)
        tmp = np.array([
            [cosa, -sina, 0.0],
            [sina, cosa, 0.0],
            [0.0, 0.0, 1.0],
        ])
        rotmat = rotmat @ tmp

    return rotmat


def _t2p_effective_distance(delta, hmin, hmax, hvert, azimuth=0.0, dip=0.0, plunge=0.0):
    rotmat = _t2p_variogram_rotmat(hmin, hmax, hvert, azimuth, dip, plunge)
    scaled = rotmat @ np.asarray(delta, dtype=float) / hmax
    return np.linalg.norm(scaled)


def _transform_vgm_locations(loc, vgm, is3d=False):
    loc = np.asarray(loc)
    if not is3d:
        _, _, _, rot, hh, _ = parse_vgm(vgm, is3d=is3d)
        return rot.apply(loc) / hh

    hmin, hmax, hvert, azimuth, dip, plunge = _parse_vgm_ranges_angles(vgm, is3d=True)
    rotmat = _t2p_variogram_rotmat(hmin, hmax, hvert, azimuth, dip, plunge)
    return (rotmat @ loc.T).T / hmax


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
    obsloc = _transform_vgm_locations(obsloc, vgm, is3d=is3d)
    newloc = _transform_vgm_locations(newloc, vgm, is3d=is3d)

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
       newloc, vgm0=_VGM_P, vgm1=_VGM_S, vgmc=_VGM_X, is3d=False,
       align_means=True, local_mean_correction=False):
    """
    2-D ordinary cokriging with a single unbiasedness constraint
    (Σ primary weights + Σ secondary weights = 1), matching the retired
    Texture2Par ordinary-cokriging reference formulation.

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

    # Apply shared rotation and range scaling. In 3D this mirrors the
    # T2P-style metric used by the 3D OK reference tests.
    obsloc0 = _transform_vgm_locations(obsloc0, vgm0, is3d=is3d)
    obsloc1 = _transform_vgm_locations(obsloc1, vgm0, is3d=is3d)
    newloc  = _transform_vgm_locations(newloc,  vgm0, is3d=is3d)

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
    # Texture2Par normally shifts secondary values to the primary global mean
    # before applying weights. Some callers pass values already aligned so a
    # target-neighbourhood reference can preserve the global, rather than local,
    # shift used by Texture2Par.
    if align_means:
        obsval1 = obsval1 + obsval0.mean() - obsval1.mean()
    est = obsval0 @ w[:n1] + obsval1 @ w[n1:n1+n2]
    if local_mean_correction:
        est += (obsval0.mean() - obsval1.mean()) * w[n1:n1+n2].sum(axis=0)
    var = np.maximum(0.0,
                     cov0_p
                     - (w[:n1+n2] * B[:n1+n2]).sum(axis=0)
                     - w[ntot-1])
    return est, var


def ck_nearest_reference(obsloc0, obsloc1, obsval0, obsval1, newloc,
                         vgm0, vgm1, vgmc, nnear,
                         local_mean_correction=False):
    """Two-variable isotropic CK using independent nearest sets per target."""
    obsloc0 = np.asarray(obsloc0)
    obsloc1 = np.asarray(obsloc1)
    newloc = np.asarray(newloc)
    estimates = np.empty(len(newloc))
    for i, target in enumerate(newloc):
        idx0 = np.argsort(np.linalg.norm(obsloc0 - target, axis=1))[:nnear]
        idx1 = np.argsort(np.linalg.norm(obsloc1 - target, axis=1))[:nnear]
        estimate, _ = ck(
            obsloc0[idx0], obsloc1[idx1], obsval0[idx0], obsval1[idx1],
            target[None, :], vgm0, vgm1, vgmc,
            align_means=False,
            local_mean_correction=local_mean_correction,
        )
        estimates[i] = estimate[0]
    return estimates


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
             options_extra="", method="OK",
             ck_mean_mode="LEGACY_STANDARDIZED") -> str:
    primary_class = classes.split()[0]
    ck_mean_line = ""
    if ck_mean_mode is not None:
        ck_mean_line = f"CK_MEAN_MODE {ck_mean_mode}"
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INTERP_METHOD {method}
          {ck_mean_line}
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_VARIANCE
          INFER_LAST_CLASS
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


def calc_ck_correlation(run_t2p, butte_ck_data, vgm_p=_VGM_P, vgm_s=_VGM_S,
                        vgm_x=_VGM_X, method="OK",
                        ck_mean_mode="LEGACY_STANDARDIZED"):
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
    t0 = time.perf_counter()
    r = run_t2p(
        _t2p_ck(vgm_s=vgm_s, vgm_x=vgm_x, vgm_p=vgm_p, method=method,
                ck_mean_mode=ck_mean_mode,
                options_extra="MAX_SECONDARY_DIST 1e10"),
        extra_files=extra,
    )
    t1 = time.perf_counter()

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
        vgm0=vgm_p,
        vgm1=vgm_s,
        vgmc=vgm_x,
    )
    t2 = time.perf_counter()

    print(f"CK: {t1 - t0:.2f}s (T2P), {t2 - t1:.2f}s (Python)")
    est = np.clip(est, 0, 1)
    corr_est = np.corrcoef(est, t2p["Layer1"].values)[0, 1]
    corr_var = np.corrcoef(var, t2p["Var_Layer1"].values)[0, 1]
    max_abs_est = np.max(np.abs(est - t2p["Layer1"].values))
    return corr_est, corr_var, max_abs_est


# ---------------------------------------------------------------------------
# 3D Reference Verification Helper
# ---------------------------------------------------------------------------


def ck3d(obsloc0, obsloc1, obsval0, obsval1,
       newloc, vgm0=_VGM_P_3D, vgm1=_VGM_S_3D, vgmc=_VGM_X_3D, is3d=True):
    return ck(obsloc0, obsloc1, obsval0, obsval1, newloc, vgm0, vgm1, vgmc, is3d)

def _override_3d_hvert(vgm, hvert_override):
    if hvert_override is None:
        return vgm
    parts = vgm.split()
    if len(parts) < 9:
        raise ValueError(f"Expected 3D variogram format, got: {vgm!r}")
    parts[5] = f"{hvert_override:g}"
    return " ".join(parts)


def _override_3d_dip(vgm, dip):
    parts = vgm.split()
    if len(parts) < 9:
        raise ValueError(f"Expected 3D variogram format, got: {vgm!r}")
    parts[7] = f"{dip:g}"
    return " ".join(parts)


def _ck_covariance_blocks(obsloc0, obsloc1, newloc,
                          vgm0=_VGM_P_3D, vgm1=_VGM_S_3D,
                          vgmc=_VGM_X_3D, is3d=True):
    obsloc0 = np.asarray(obsloc0, dtype=float)
    obsloc1 = np.asarray(obsloc1, dtype=float)
    newloc = np.asarray(newloc, dtype=float)

    vtype0, nugget0, sill0, rot0, hh0, _ = parse_vgm(vgm0, is3d=is3d)
    vtype1, nugget1, sill1, rot1, hh1, _ = parse_vgm(vgm1, is3d=is3d)
    vtypec, _, sillc, rotc, hhc, _ = parse_vgm(vgmc, is3d=is3d)
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
        newloc = np.column_stack([newloc, np.zeros(len(newloc))])

    loc0 = _transform_vgm_locations(obsloc0, vgm0, is3d=is3d)
    loc1 = _transform_vgm_locations(obsloc1, vgm0, is3d=is3d)
    targets = _transform_vgm_locations(newloc, vgm0, is3d=is3d)

    C11 = _covfunc[vtype0](squareform(pdist(loc0))) * sill0
    np.fill_diagonal(C11, nugget0 + sill0)
    C22 = _covfunc[vtype1](squareform(pdist(loc1))) * sill1
    np.fill_diagonal(C22, nugget1 + sill1)
    C12 = _covfunc[vtypec](cdist(loc0, loc1)) * sillc
    c10 = _covfunc[vtype0](cdist(loc0, targets)) * sill0
    c20 = _covfunc[vtypec](cdist(loc1, targets)) * sillc
    return {"C11": C11, "C22": C22, "C12": C12, "c10": c10, "c20": c20}


def _ck_block_covariance_summaries(obsloc0, obsloc1, newloc,
                                   vgm0=_VGM_P_3D, vgm1=_VGM_S_3D,
                                   vgmc=_VGM_X_3D, dip0=0.0, dip1=5.0):
    blocks0 = _ck_covariance_blocks(
        obsloc0,
        obsloc1,
        newloc,
        _override_3d_dip(vgm0, dip0),
        _override_3d_dip(vgm1, dip0),
        _override_3d_dip(vgmc, dip0),
        is3d=True,
    )
    blocks1 = _ck_covariance_blocks(
        obsloc0,
        obsloc1,
        newloc,
        _override_3d_dip(vgm0, dip1),
        _override_3d_dip(vgm1, dip1),
        _override_3d_dip(vgmc, dip1),
        is3d=True,
    )
    summaries = {}
    for name in ("C11", "C22", "C12", "c10", "c20"):
        diff = np.abs(blocks1[name] - blocks0[name])
        summaries[name] = {
            "dip0_min": float(blocks0[name].min()),
            "dip0_mean": float(blocks0[name].mean()),
            "dip0_max": float(blocks0[name].max()),
            "dip1_min": float(blocks1[name].min()),
            "dip1_mean": float(blocks1[name].mean()),
            "dip1_max": float(blocks1[name].max()),
            "mean_abs_delta": float(diff.mean()),
            "max_abs_delta": float(diff.max()),
        }
    return summaries


def _ck_rhs_covariance_diagnostics(obsloc0, obsloc1, newloc,
                                   vgm0=_VGM_P_3D, vgmc=_VGM_X_3D,
                                   is3d=True):
    obsloc0 = np.asarray(obsloc0)
    obsloc1 = np.asarray(obsloc1)
    newloc = np.asarray(newloc)

    vtype0, _, sill0, rot0, hh0, _ = parse_vgm(vgm0, is3d=is3d)
    vtypec, _, sillc, rotc, hhc, _ = parse_vgm(vgmc, is3d=is3d)
    assert np.allclose(rot0.as_matrix(), rotc.as_matrix(), atol=1e-6), \
        "All variograms must share the same rotation"
    assert np.allclose(hh0, hhc, atol=1e-6), \
        "All variograms must share the same range parameters"

    if not is3d:
        obsloc0 = np.column_stack([obsloc0, np.zeros(len(obsloc0))])
        obsloc1 = np.column_stack([obsloc1, np.zeros(len(obsloc1))])
        newloc = np.column_stack([newloc, np.zeros(len(newloc))])

    obsloc0 = _transform_vgm_locations(obsloc0, vgm0, is3d=is3d)
    obsloc1 = _transform_vgm_locations(obsloc1, vgm0, is3d=is3d)
    newloc = _transform_vgm_locations(newloc, vgm0, is3d=is3d)

    c10 = _covfunc[vtype0](cdist(obsloc0, newloc)) * sill0
    c20 = _covfunc[vtypec](cdist(obsloc1, newloc)) * sillc
    return c10.max(axis=0), c20.max(axis=0)


def _print_series_summary(label, values):
    values = np.asarray(values, dtype=float)
    print(
        f"{label} min/mean/max: "
        f"{values.min():.8g} / {values.mean():.8g} / {values.max():.8g}"
    )


def _print_comparison_diagnostics(label, targets, ref, t2p, include_rows=True):
    ref = np.asarray(ref, dtype=float)
    t2p = np.asarray(t2p, dtype=float)
    diff = t2p - ref
    rmse = float(np.sqrt(np.mean(diff * diff)))
    if np.std(ref) == 0.0 or np.std(t2p) == 0.0:
        corr = np.nan
    else:
        corr = float(np.corrcoef(ref, t2p)[0, 1])

    print(f"{label} reference min/max: {ref.min():.8g} / {ref.max():.8g}")
    print(f"{label} T2P min/max:       {t2p.min():.8g} / {t2p.max():.8g}")
    print(f"{label} max abs diff:     {np.max(np.abs(diff)):.8g}")
    print(f"{label} RMSE:             {rmse:.8g}")
    print(f"{label} correlation:      {corr:.8g}")

    if include_rows:
        diag = targets[["X", "Y", "Z"]].copy()
        diag["ref_est" if label == "estimate" else "ref_var"] = ref
        diag["t2p_est" if label == "estimate" else "t2p_var"] = t2p
        diag["diff"] = diff
        worst = diag.iloc[np.argsort(np.abs(diff))[::-1][:10]]
        print(f"{label} worst 10 rows:")
        print(worst.to_string(index=False))


def _comparison_stats(ref, t2p):
    ref = np.asarray(ref, dtype=float)
    t2p = np.asarray(t2p, dtype=float)
    diff = t2p - ref
    corr = np.nan
    if np.std(ref) > 0.0 and np.std(t2p) > 0.0:
        corr = float(np.corrcoef(ref, t2p)[0, 1])
    return {
        "max_abs_diff": float(np.max(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff * diff))),
        "corr": corr,
        "ref_min": float(ref.min()),
        "ref_max": float(ref.max()),
        "t2p_min": float(t2p.min()),
        "t2p_max": float(t2p.max()),
    }


def _print_stats(label, stats):
    print(
        f"{label}: max_abs_diff={stats['max_abs_diff']:.8g}, "
        f"rmse={stats['rmse']:.8g}, corr={stats['corr']:.8g}, "
        f"ref_min/max={stats['ref_min']:.8g}/{stats['ref_max']:.8g}, "
        f"t2p_min/max={stats['t2p_min']:.8g}/{stats['t2p_max']:.8g}"
    )


# The Python CK reference uses the full 150-point secondary subsample for
# every target. T2P has a finite default MAX_SECONDARY_DIST, so these
# full-matrix reference tests explicitly set MAX_SECONDARY_DIST large enough
# to prevent distance-based secondary filtering.
def _t2p_ck3d(classes="Coarse\n  Fine", covariate="Rho",
             datafile_p="pc3d.csv", datafile_s="aem3d.csv",
             gridfile="grid.grd",
             vgm_p=_VGM_P_3D, vgm_s=_VGM_S_3D, vgm_x=_VGM_X_3D,
             transform_s="LOG10",
             options_extra=None):
    if options_extra is None:
        options_extra = """\
          INTERP_DIM          3
          WRITE_DATASET_FILES
          MAX_SECONDARY_DIST  1e10
        """
    return _t2p_ck(
        classes, covariate, datafile_p, datafile_s, gridfile,
        vgm_p, vgm_s, vgm_x, transform_s, options_extra,
        method="OK", ck_mean_mode="LEGACY_STANDARDIZED",
    )


CK3D_SYNTHETIC_GRID = """\
NX    2
NY    2
NZ    2
DX    10.0
DY    10.0
TOP   constant 100.0
BOT   constant 80.0
IBOUND constant 1
"""

CK3D_SYNTHETIC_PRIMARY = """\
Well,ID,n,x,y,Zland,Depth0,Depth1,Coarse
P1,1,1,2.0,2.0,100.0,4.0,6.0,0.22
P2,2,1,18.0,3.0,100.0,8.0,10.0,0.38
P3,3,1,4.0,17.0,100.0,13.0,15.0,0.72
P4,4,1,17.0,18.0,100.0,17.0,19.0,0.84
P5,5,1,9.0,11.0,100.0,9.0,11.0,0.55
"""

CK3D_SYNTHETIC_SECONDARY = """\
iRho,ID,n,x,y,Zland,Depth0,Depth1,Rho
S1,1,1,1.0,12.0,100.0,5.0,7.0,8.0
S2,2,1,7.0,4.0,100.0,11.0,13.0,11.0
S3,3,1,14.0,8.0,100.0,7.0,9.0,16.0
S4,4,1,19.0,14.0,100.0,15.0,17.0,23.0
S5,5,1,11.0,19.0,100.0,3.0,5.0,31.0
S6,6,1,5.0,18.0,100.0,18.0,20.0,19.0
"""


def _shift_xy_dataset(text, dx, dy):
    df = pd.read_csv(StringIO(text))
    x_col = "x" if "x" in df.columns else "X"
    y_col = "y" if "y" in df.columns else "Y"
    df[x_col] = df[x_col] + dx
    df[y_col] = df[y_col] + dy
    return df.to_csv(index=False, lineterminator="\n")


def _t2p_ck3d_synthetic(vgm_p, vgm_s, vgm_x, xoffset=None, yoffset=None):
    offset_lines = ""
    if xoffset is not None:
        offset_lines += f"          XOFFSET   {xoffset}\n"
    if yoffset is not None:
        offset_lines += f"          YOFFSET   {yoffset}\n"
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INTERP_METHOD OK
          CK_MEAN_MODE LEGACY_STANDARDIZED
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_VARIANCE
          INFER_LAST_CLASS
          INTERP_DIM 3
          WRITE_DATASET_FILES
          MAX_SECONDARY_DIST 1e10
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN COVARIATES
          Coarse Rho
        END COVARIATES

        BEGIN DATASET
          FILE pc3d.csv
        END DATASET

        BEGIN DATASET
          FILE aem3d.csv
          TRANSFORM LOG10
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
{offset_lines.rstrip()}
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 {vgm_p}
          CLASS Rho
            1 {vgm_s}
          CLASS Coarse:Rho
            1 {vgm_x}
        END VARIOGRAMS
        """


def _synthetic_ck3d_reference_inputs(primary_text=CK3D_SYNTHETIC_PRIMARY,
                                     secondary_text=CK3D_SYNTHETIC_SECONDARY):
    primary = pd.read_csv(StringIO(primary_text))
    secondary = pd.read_csv(StringIO(secondary_text))
    primary["z"] = primary["Zland"] - 0.5 * (primary["Depth0"] + primary["Depth1"])
    secondary["z"] = secondary["Zland"] - 0.5 * (secondary["Depth0"] + secondary["Depth1"])

    s_log10 = np.log10(secondary["Rho"].values)
    v_adjust = primary["Coarse"].mean() - s_log10.mean()
    s_val = s_log10 + v_adjust
    return primary, secondary, s_val


def _synthetic_ck3d_output_arrays(output, layer_z=(95.0, 85.0)):
    targets = []
    est = []
    var = []
    for row in output.itertuples(index=False):
        for layer_index, z in enumerate(layer_z, start=1):
            targets.append([row.X, row.Y, z])
            est.append(getattr(row, f"Layer{layer_index}"))
            var.append(getattr(row, f"Var_Layer{layer_index}"))
    return np.asarray(targets), np.asarray(est), np.asarray(var)


def _run_synthetic_ck3d(run_t2p, vgm_p, vgm_s, vgm_x,
                        primary_text=CK3D_SYNTHETIC_PRIMARY,
                        secondary_text=CK3D_SYNTHETIC_SECONDARY,
                        xoffset=None, yoffset=None):
    r = run_t2p(
        _t2p_ck3d_synthetic(vgm_p, vgm_s, vgm_x, xoffset=xoffset, yoffset=yoffset),
        extra_files={
            "pc3d.csv": primary_text,
            "aem3d.csv": secondary_text,
            "grid.grd": CK3D_SYNTHETIC_GRID,
        },
    )
    assert r.returncode == 0, r.combined
    t2p = pd.read_csv(r.workdir / "test_Coarse.csv")
    targets, t2p_est, t2p_var = _synthetic_ck3d_output_arrays(t2p)
    primary, secondary, s_val = _synthetic_ck3d_reference_inputs(primary_text, secondary_text)
    ref_est, ref_var = ck3d(
        primary[["x", "y", "z"]].values,
        secondary[["x", "y", "z"]].values,
        primary["Coarse"].values,
        s_val,
        targets,
        vgm0=vgm_p,
        vgm1=vgm_s,
        vgmc=vgm_x,
    )
    return {
        "targets": targets,
        "ref_est": np.clip(ref_est, 0, 1),
        "ref_var": ref_var,
        "t2p_est": t2p_est,
        "t2p_var": t2p_var,
        "t2p": t2p,
        "primary": primary,
        "secondary": secondary,
    }


def _calc_ck_correlation_3d_result(
    run_t2p,
    butte_ck_data,
    vgm_p=_VGM_P_3D,
    vgm_s=_VGM_S_3D,
    vgm_x=_VGM_X_3D,
    z_mode="random",
    hvert_override=None,
    options_extra=None,
    debug=False,
):
    """
    Run one 3D CK reference comparison and return both Python and T2P arrays.

    The Python reference includes all primary and secondary observations
    passed to ck3d. It does not implement T2P's target-specific
    secondary-distance filtering, so exact reference comparisons must run T2P
    with MAX_SECONDARY_DIST large enough to include the selected secondary
    subsample for every target.
    """
    vgm_p = _override_3d_hvert(vgm_p, hvert_override)
    vgm_s = _override_3d_hvert(vgm_s, hvert_override)
    vgm_x = _override_3d_hvert(vgm_x, hvert_override)

    obs = pd.read_csv(StringIO(butte_ck_data["files"]["pc2d.csv"]))
    aem_all = pd.read_csv(StringIO(butte_ck_data["files"]["aem2d.csv"]))

    rng = np.random.default_rng(42)
    idx = rng.choice(len(aem_all), 150, replace=False)
    aem_sub = aem_all.iloc[idx].reset_index(drop=True)
    aem_sub["ID"] = np.arange(1, len(aem_sub) + 1)

    if z_mode == "random":
        dz = 75
        obs["Zland"] = rng.uniform(-dz, dz, len(obs)) + 52
        aem_sub["Zland"] = rng.uniform(-dz, dz, len(aem_sub)) + 52
    elif z_mode == "flat50":
        obs["Zland"] = 52.0
        aem_sub["Zland"] = 52.0
    else:
        raise ValueError(f"Unknown z_mode {z_mode!r}")
    obs["Depth0"] = 1.0
    obs["Depth1"] = 3.0
    aem_sub["Depth0"] = 1.0
    aem_sub["Depth1"] = 3.0

    extra = {
        **butte_ck_data["files"],
        "pc3d.csv": obs.to_csv(index=False, lineterminator='\n'),
        "aem3d.csv": aem_sub.to_csv(index=False, lineterminator='\n')
    }

    t0 = time.perf_counter()
    r = run_t2p(
        _t2p_ck3d(
            vgm_s=vgm_s,
            vgm_x=vgm_x,
            vgm_p=vgm_p,
            options_extra=options_extra,
        ),
        extra_files=extra,
    )
    t1 = time.perf_counter()
    assert r.returncode == 0, r.combined

    t2p = pd.read_csv(r.workdir / "test_Coarse.csv")
    t2p["Z"] = 50.0

    # T2P shifts secondary values so they share the primary mean before
    # multiplying by weights (Interpolators.f90, line ~221):
    #   v_adjust = mean(primary) - mean(secondary)
    # Apply the same adjustment here so the Python estimate matches T2P.
    s_log10 = np.log10(aem_sub["Rho"].values)
    v_adjust = obs["Coarse"].mean() - s_log10.mean()
    s_val = s_log10 + v_adjust

    # Perform 3D validation step
    obs    ["z"] = obs    ["Zland"] - 2
    aem_sub["z"] = aem_sub["Zland"] - 2
    est, var = ck3d(
        obs[["x", "y", "z"]].values,
        aem_sub[["x", "y", "z"]].values,
        obs["Coarse"].values,
        s_val,
        t2p[["X", "Y", "Z"]].values,
        vgm0=vgm_p,
        vgm1=vgm_s,
        vgmc=vgm_x,
    )
    t2 = time.perf_counter()

    est = np.clip(est, 0, 1)
    if debug:
        print(f"CK 3D: {t1 - t0:.2f}s (T2P), {t2 - t1:.2f}s (Python)")
        print(f"CK 3D diagnostic z_mode={z_mode}, hvert_override={hvert_override}")
        print(f"vgm_p: {vgm_p}")
        print(f"vgm_s: {vgm_s}")
        print(f"vgm_x: {vgm_x}")
        _print_comparison_diagnostics(
            "estimate",
            t2p[["X", "Y", "Z"]],
            est,
            t2p["Layer1"].values,
        )
        _print_comparison_diagnostics(
            "variance",
            t2p[["X", "Y", "Z"]],
            var,
            t2p["Var_Layer1"].values,
        )
        c10_max, c20_max = _ck_rhs_covariance_diagnostics(
            obs[["x", "y", "z"]].values,
            aem_sub[["x", "y", "z"]].values,
            t2p[["X", "Y", "Z"]].values,
            vgm0=vgm_p,
            vgmc=vgm_x,
            is3d=True,
        )
        _print_series_summary("c10 max by target", c10_max)
        _print_series_summary("c20 max by target", c20_max)

    return {
        "est": est,
        "var": var,
        "t2p_est": t2p["Layer1"].values,
        "t2p_var": t2p["Var_Layer1"].values,
        "t2p": t2p,
        "obs": obs,
        "aem_sub": aem_sub,
        "vgm_p": vgm_p,
        "vgm_s": vgm_s,
        "vgm_x": vgm_x,
    }


def calc_ck_correlation_3d(
    run_t2p,
    butte_ck_data,
    vgm_p=_VGM_P_3D,
    vgm_s=_VGM_S_3D,
    vgm_x=_VGM_X_3D,
    z_mode="random",
    hvert_override=None,
    debug=False,
):
    """
    3D CK reference calculation generating random Z coordinates to test cross-validation.
    """
    result = _calc_ck_correlation_3d_result(
        run_t2p,
        butte_ck_data,
        vgm_p=vgm_p,
        vgm_s=vgm_s,
        vgm_x=vgm_x,
        z_mode=z_mode,
        hvert_override=hvert_override,
        debug=debug,
    )
    est = result["est"]
    var = result["var"]
    t2p_est = result["t2p_est"]
    t2p_var = result["t2p_var"]

    np.testing.assert_allclose(est, t2p_est, atol=1e-5)
    np.testing.assert_allclose(var, t2p_var, atol=1e-5)
    corr_est = np.corrcoef(est, t2p_est)[0, 1]
    corr_var = np.corrcoef(var, t2p_var)[0, 1]
    assert corr_est > 0.99, (
        f"CK estimate correlation with reference = {corr_est:.4f} (expected > 0.99)"
    )
    assert corr_var > 0.99, (
        f"CK variance correlation with reference = {corr_var:.4f} (expected > 0.99)"
    )


# ---------------------------------------------------------------------------
# 1. Cokriging tests (primary + secondary AEM dataset)
# ---------------------------------------------------------------------------

class TestCokrigingButteCK:

    def test_ck_nthreads_12_matches_single_thread(self, run_t2p, butte_ck_data):
        single = run_t2p(
            _t2p_ck(options_extra="NTHREADS 1"),
            filename="single.t2p",
            extra_files=butte_ck_data["files"],
        )
        assert single.returncode == 0, single.combined
        single_output = pd.read_csv(single.workdir / "single_Coarse.csv").copy()

        parallel = run_t2p(
            _t2p_ck(options_extra="NTHREADS 12"),
            filename="parallel.t2p",
            extra_files=butte_ck_data["files"],
        )
        assert parallel.returncode == 0, parallel.combined
        parallel_output = pd.read_csv(parallel.workdir / "parallel_Coarse.csv")

        assert parallel_output["Layer1"].to_numpy() == pytest.approx(
            single_output["Layer1"].to_numpy(), abs=1e-12
        )
        assert parallel_output["Var_Layer1"].to_numpy() == pytest.approx(
            single_output["Var_Layer1"].to_numpy(), abs=1e-12
        )

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

    def test_legacy_standardized_ck_reference_correlation(self, run_t2p, butte_ck_data):
        """
        Libkriging CK with historical Texture2Par mean standardization matches
        the isotropic full-secondary Python reference.
        """

        corr_est, corr_var, max_abs_est = calc_ck_correlation(run_t2p, butte_ck_data)
        assert corr_est > 0.99, (
            f"CK estimate correlation with reference = {corr_est:.4f} (expected > 0.99)"
        )
        assert corr_var > 0.99, (
            f"CK variance correlation with reference = {corr_var:.4f} (expected > 0.99)"
        )
        assert max_abs_est < 1e-5, (
            f"CK maximum estimate difference = {max_abs_est:.8g} (expected < 1e-5)"
        )

    def test_ck_mean_modes_match_references_and_differ(
        self, run_t2p, butte_ck_data
    ):
        nnear = 5
        vgm_p = "sph 0.0 0.12 5000 5000 0.0 5"
        vgm_s = "sph 0.0 0.12 5000 5000 0.0 5"
        vgm_x = "sph 0.0 0.048 5000 5000 0.0"
        default = run_t2p(
            _t2p_ck(method="OK", vgm_p=vgm_p, vgm_s=vgm_s, vgm_x=vgm_x,
                    ck_mean_mode=None,
                    options_extra="MAX_SECONDARY_DIST 1e10"),
            extra_files=butte_ck_data["files"],
        )
        assert default.returncode == 0, default.combined
        default_output = pd.read_csv(
            default.workdir / "test_Coarse.csv"
        ).head(25).copy()
        default_mode_line = next(
            line for line in default.log.splitlines() if "CK Mean Mode:" in line
        )
        assert default_mode_line.rstrip().endswith("LIBKRIGING")

        legacy = run_t2p(
            _t2p_ck(method="OK", vgm_p=vgm_p, vgm_s=vgm_s, vgm_x=vgm_x,
                    ck_mean_mode="LEGACY_STANDARDIZED",
                    options_extra="MAX_SECONDARY_DIST 1e10"),
            extra_files=butte_ck_data["files"],
        )
        assert legacy.returncode == 0, legacy.combined
        legacy_output = pd.read_csv(
            legacy.workdir / "test_Coarse.csv"
        ).head(25).copy()
        legacy_mode_line = next(
            line for line in legacy.log.splitlines() if "CK Mean Mode:" in line
        )
        assert legacy_mode_line.rstrip().endswith("LEGACY_STANDARDIZED")

        primary = pd.read_csv(StringIO(butte_ck_data["files"]["pc2d.csv"]))
        secondary = pd.read_csv(StringIO(butte_ck_data["files"]["aem2d.csv"]))
        primary_values = primary["Coarse"].to_numpy()
        raw_secondary_values = np.log10(secondary["Rho"].to_numpy())
        legacy_secondary_values = raw_secondary_values.copy()
        legacy_secondary_values += (
            primary_values.mean() - legacy_secondary_values.mean()
        )

        default_reference = ck_nearest_reference(
            primary[["x", "y"]].to_numpy(),
            secondary[["x", "y"]].to_numpy(),
            primary_values,
            raw_secondary_values,
            default_output[["X", "Y"]].to_numpy(),
            vgm_p, vgm_s, vgm_x, nnear,
            local_mean_correction=True,
        )
        legacy_reference = ck_nearest_reference(
            primary[["x", "y"]].to_numpy(),
            secondary[["x", "y"]].to_numpy(),
            primary_values,
            legacy_secondary_values,
            legacy_output[["X", "Y"]].to_numpy(),
            vgm_p, vgm_s, vgm_x, nnear,
        )
        default_actual = default_output["Layer1"].to_numpy()
        legacy_actual = legacy_output["Layer1"].to_numpy()
        assert default_actual == pytest.approx(
            np.clip(default_reference, 0.0, 1.0), abs=1e-5
        )
        assert legacy_actual == pytest.approx(
            np.clip(legacy_reference, 0.0, 1.0), abs=1e-5
        )
        assert np.max(np.abs(default_actual - legacy_actual)) > 1e-3

    def test_libkriging_ck_respects_max_secondary_distance(
        self, run_t2p, butte_ck_data
    ):
        full = run_t2p(
            _t2p_ck(method="OK", ck_mean_mode="LIBKRIGING",
                    options_extra="MAX_SECONDARY_DIST 1e10"),
            extra_files=butte_ck_data["files"],
        )
        assert full.returncode == 0, full.combined
        full_values = pd.read_csv(full.workdir / "test_Coarse.csv")["Layer1"].values.copy()
        restricted = run_t2p(
            _t2p_ck(method="OK", ck_mean_mode="LIBKRIGING",
                    options_extra="MAX_SECONDARY_DIST 100"),
            extra_files=butte_ck_data["files"],
        )
        assert restricted.returncode == 0, restricted.combined
        restricted_values = pd.read_csv(
            restricted.workdir / "test_Coarse.csv"
        )["Layer1"].values
        assert np.max(np.abs(full_values - restricted_values)) > 1e-5

    def test_libkriging_ck_supports_multiple_secondary_variables(
        self, run_t2p, butte_ck_data
    ):
        rho = pd.read_csv(StringIO(butte_ck_data["files"]["aem2d.csv"])).head(8)
        sal = rho.drop(columns=["Rho"]).copy()
        sal["x"] += 10.0
        sal["Sal"] = np.linspace(0.2, 0.9, len(sal))
        t2p = _t2p_ck(method="OK_LIB", ck_mean_mode="LIBKRIGING",
                      options_extra="MAX_SECONDARY_DIST 1e10")
        t2p = t2p.replace(
            "Coarse Rho\n        END COVARIATES",
            "Coarse Rho Sal\n        END COVARIATES",
        )
        t2p = t2p.replace(
            "        BEGIN FLOW_MODEL",
            "        BEGIN DATASET\n"
            "          FILE sal2d.csv\n"
            "        END DATASET\n\n"
            "        BEGIN FLOW_MODEL",
            1,
        )
        t2p = t2p.replace(
            "        END VARIOGRAMS",
            "          CLASS Sal\n"
            "            1 sph 0.0 0.20 5000 5000 0.0 8\n"
            "          CLASS Coarse:Sal\n"
            "            1 sph 0.0 0.02 5000 5000 0.0\n"
            "          CLASS Rho:Sal\n"
            "            1 sph 0.0 0.01 5000 5000 0.0\n"
            "        END VARIOGRAMS",
            1,
        )
        files = {
            **butte_ck_data["files"],
            "aem2d.csv": rho.to_csv(index=False, lineterminator="\n"),
            "sal2d.csv": sal.to_csv(index=False, lineterminator="\n"),
        }
        r = run_t2p(t2p, extra_files=files)
        assert r.returncode == 0, r.combined
        values = pd.read_csv(r.workdir / "test_Coarse.csv")["Layer1"]
        assert np.isfinite(values).all()

    @pytest.mark.parametrize("aniso1",  [1.0, 0.1])
    @pytest.mark.parametrize("azimuth", [0.0, 225])
    def test_ck_reference_correlation_aniso(self, run_t2p, butte_ck_data, aniso1, azimuth):
        """
        CK estimates and variances correlate > 0.99 with a Python reference with anisotropic variograms.
        """
        vgm_p=f"sph 0.0    0.12  {5000*aniso1} 5000 {azimuth} 150"
        vgm_s=f"sph 0.02   0.08  {5000*aniso1} 5000 {azimuth} 150"
        vgm_x=f"sph 0.0   0.048  {5000*aniso1} 5000 {azimuth}"
        corr_est, corr_var, _ = calc_ck_correlation(
            run_t2p, butte_ck_data, vgm_p, vgm_s, vgm_x,
            method="OK", ck_mean_mode="LEGACY_STANDARDIZED",
        )
        assert corr_est > 0.99, (
            f"CK estimate correlation with reference = {corr_est:.4f} (expected > 0.99)"
        )
        assert corr_var > 0.99, (
            f"CK variance correlation with reference = {corr_var:.4f} (expected > 0.99)"
        )



class TestCokrigingButteCK3D:

    def _zero_rotation_vgms(self):
        return (
            "sph 0.0    0.12  500 5000 500 0.0 0.0 0.0 150",
            "sph 0.02   0.08  500 5000 500 0.0 0.0 0.0 150",
            "sph 0.0   0.048  500 5000 500 0.0 0.0 0.0",
        )

    def _synthetic_dipped_vgms(self):
        return (
            "sph 0.0    0.12  12.0 25.0 18.0 0.0 5.0 0.0 20",
            "sph 0.02   0.08  12.0 25.0 18.0 0.0 5.0 0.0 20",
            "sph 0.0    0.04  12.0 25.0 18.0 0.0 5.0 0.0",
        )

    def test_ck3d_rotation_metric_uses_t2p_fortran_dip_convention(self):
        delta = np.array([-1.0, 15.0, -9.0])
        hmin, hmax, hvert = 500.0, 5000.0, 500.0

        zero = _t2p_variogram_rotmat(hmin, hmax, hvert, 0.0, 0.0, 0.0)
        np.testing.assert_allclose(
            zero,
            np.diag([hmax / hmin, 1.0, hmax / hvert]),
            atol=1e-14,
        )

        azimuth = _t2p_variogram_rotmat(hmin, hmax, hvert, 35.0, 0.0, 0.0)
        assert _t2p_effective_distance(delta, hmin, hmax, hvert, 35.0, 0.0, 0.0) == pytest.approx(
            np.linalg.norm((azimuth @ delta) / hmax),
            abs=1e-14,
        )

        dip = _t2p_variogram_rotmat(12.0, 25.0, 18.0, 0.0, 5.0, 0.0)
        expected = np.array([
            [25.0 / 12.0, 0.0, 0.0],
            [0.0, np.cos(np.deg2rad(5.0)), -np.sin(np.deg2rad(5.0))],
            [0.0, (25.0 / 18.0) * np.sin(np.deg2rad(5.0)), (25.0 / 18.0) * np.cos(np.deg2rad(5.0))],
        ])
        np.testing.assert_allclose(dip, expected, atol=1e-14)

    def test_ck3d_dip_pair_trace_matches_fortran_covariance(self):
        hmin, hmax, hvert = 12.0, 25.0, 18.0
        rotmat = _t2p_variogram_rotmat(
            hmin, hmax, hvert,
            azimuth=0.0,
            dip=5.0,
            plunge=0.0,
        )

        cases = [
            (np.array([-1.0, 15.0, -9.0]), 0.12, 0.00923520489567),
            (np.array([-4.0, 4.0, 5.0]), 0.12, 0.0419239732925),
        ]
        for delta, sill, expected_cov in cases:
            dist = np.linalg.norm(rotmat @ delta)
            hr = dist / hmax
            cov = sill * (1.0 - 1.5 * hr + 0.5 * hr**3) if dist < hmax else 0.0
            assert cov == pytest.approx(expected_cov, abs=1e-12)

    def test_ck3d_t2p_style_metric_is_common_translation_invariant(self):
        vgm_p, vgm_s, vgm_x = self._synthetic_dipped_vgms()
        primary, secondary, _ = _synthetic_ck3d_reference_inputs()
        targets = np.array([
            [5.0, 5.0, 95.0],
            [15.0, 5.0, 95.0],
            [5.0, 15.0, 85.0],
            [15.0, 15.0, 85.0],
        ])
        shift = np.array([10000.0, -20000.0, 123.0])

        blocks0 = _ck_covariance_blocks(
            primary[["x", "y", "z"]].values,
            secondary[["x", "y", "z"]].values,
            targets,
            vgm_p,
            vgm_s,
            vgm_x,
            is3d=True,
        )
        blocks1 = _ck_covariance_blocks(
            primary[["x", "y", "z"]].values + shift,
            secondary[["x", "y", "z"]].values + shift,
            targets + shift,
            vgm_p,
            vgm_s,
            vgm_x,
            is3d=True,
        )

        for name in ("C11", "C22", "C12", "c10", "c20"):
            np.testing.assert_allclose(blocks0[name], blocks1[name], atol=1e-14)

    def test_ck3d_python_block_covariance_summaries_are_dip_sensitive(self):
        vgm_p, vgm_s, vgm_x = self._synthetic_dipped_vgms()
        primary, secondary, _ = _synthetic_ck3d_reference_inputs()
        targets = np.array([
            [5.0, 5.0, 95.0],
            [15.0, 5.0, 95.0],
            [5.0, 15.0, 85.0],
            [15.0, 15.0, 85.0],
        ])

        summaries = _ck_block_covariance_summaries(
            primary[["x", "y", "z"]].values,
            secondary[["x", "y", "z"]].values,
            targets,
            vgm_p,
            vgm_s,
            vgm_x,
        )
        assert max(s["max_abs_delta"] for s in summaries.values()) > 1e-4
        assert summaries["c20"]["max_abs_delta"] > 1e-5

    def test_ck3d_dip_synthetic_local_coordinates_matches_reference(self, run_t2p):
        vgm_p, vgm_s, vgm_x = self._synthetic_dipped_vgms()
        result = _run_synthetic_ck3d(run_t2p, vgm_p, vgm_s, vgm_x)

        np.testing.assert_allclose(result["ref_est"], result["t2p_est"], atol=1e-5)
        np.testing.assert_allclose(result["ref_var"], result["t2p_var"], atol=1e-5)

    def test_ck3d_dip_synthetic_translation_invariant(self, run_t2p):
        vgm_p, vgm_s, vgm_x = self._synthetic_dipped_vgms()
        local = _run_synthetic_ck3d(run_t2p, vgm_p, vgm_s, vgm_x)
        dx, dy = 10000.0, 20000.0
        shifted = _run_synthetic_ck3d(
            run_t2p,
            vgm_p,
            vgm_s,
            vgm_x,
            primary_text=_shift_xy_dataset(CK3D_SYNTHETIC_PRIMARY, dx, dy),
            secondary_text=_shift_xy_dataset(CK3D_SYNTHETIC_SECONDARY, dx, dy),
            xoffset=dx,
            yoffset=dy,
        )

        np.testing.assert_allclose(local["t2p_est"], shifted["t2p_est"], atol=1e-5)
        np.testing.assert_allclose(local["t2p_var"], shifted["t2p_var"], atol=1e-5)

    def test_ck3d_flat_z_matches_reference_zero_rotation(self, run_t2p, butte_ck_data):
        vgm_p, vgm_s, vgm_x = self._zero_rotation_vgms()

        calc_ck_correlation_3d(
            run_t2p,
            butte_ck_data,
            vgm_p=vgm_p,
            vgm_s=vgm_s,
            vgm_x=vgm_x,
            z_mode="flat50",
        )

    def test_ck3d_random_z_huge_hvert_matches_reference_zero_rotation(self, run_t2p, butte_ck_data):
        vgm_p, vgm_s, vgm_x = self._zero_rotation_vgms()

        calc_ck_correlation_3d(
            run_t2p,
            butte_ck_data,
            vgm_p=vgm_p,
            vgm_s=vgm_s,
            vgm_x=vgm_x,
            z_mode="random",
            hvert_override=1e9,
        )

    def test_ck3d_random_z_finite_hvert_diagnostic_zero_rotation(self, run_t2p, butte_ck_data):
        vgm_p, vgm_s, vgm_x = self._zero_rotation_vgms()

        calc_ck_correlation_3d(
            run_t2p,
            butte_ck_data,
            vgm_p=vgm_p,
            vgm_s=vgm_s,
            vgm_x=vgm_x,
            z_mode="random",
        )

    def test_ck3d_full_secondary_reference_requires_large_max_secondary_dist(
        self,
        run_t2p,
        butte_ck_data,
    ):
        vgm_p, vgm_s, vgm_x = self._zero_rotation_vgms()

        full = _calc_ck_correlation_3d_result(
            run_t2p,
            butte_ck_data,
            vgm_p=vgm_p,
            vgm_s=vgm_s,
            vgm_x=vgm_x,
            z_mode="flat50",
        )
        np.testing.assert_allclose(full["est"], full["t2p_est"], atol=1e-5)

        # This intentionally restrictive run demonstrates the test contract:
        # the Python reference still uses the full secondary subsample, while
        # T2P solves a different, distance-filtered system.
        restricted = _calc_ck_correlation_3d_result(
            run_t2p,
            butte_ck_data,
            vgm_p=vgm_p,
            vgm_s=vgm_s,
            vgm_x=vgm_x,
            z_mode="flat50",
            options_extra="""\
              INTERP_DIM          3
              WRITE_DATASET_FILES
              MAX_SECONDARY_DIST  1.0
            """,
        )

        assert np.max(np.abs(restricted["t2p_est"] - full["t2p_est"])) > 1e-3
        assert np.max(np.abs(restricted["t2p_est"] - restricted["est"])) > 1e-3

    def test_ck3d_butte_dip_vs_no_dip_diagnostic(self, run_t2p, butte_ck_data):
        base = {
            "aniso1": 0.1,
            "aniso2": 0.1,
            "azimuth": 225.0,
        }

        def vgms(dip):
            return (
                f"sph 0.0    0.12  {5000*base['aniso1']} 5000 {5000*base['aniso2']} {base['azimuth']} {dip} 0.0 150",
                f"sph 0.02   0.08  {5000*base['aniso1']} 5000 {5000*base['aniso2']} {base['azimuth']} {dip} 0.0 150",
                f"sph 0.0   0.048  {5000*base['aniso1']} 5000 {5000*base['aniso2']} {base['azimuth']} {dip} 0.0",
            )

        vgm_p0, vgm_s0, vgm_x0 = vgms(0.0)
        no_dip = _calc_ck_correlation_3d_result(
            run_t2p,
            butte_ck_data,
            vgm_p=vgm_p0,
            vgm_s=vgm_s0,
            vgm_x=vgm_x0,
            z_mode="random",
        )
        np.testing.assert_allclose(no_dip["est"], no_dip["t2p_est"], atol=1e-5)
        np.testing.assert_allclose(no_dip["var"], no_dip["t2p_var"], atol=1e-5)

        vgm_p5, vgm_s5, vgm_x5 = vgms(5.0)
        dip = _calc_ck_correlation_3d_result(
            run_t2p,
            butte_ck_data,
            vgm_p=vgm_p5,
            vgm_s=vgm_s5,
            vgm_x=vgm_x5,
            z_mode="random",
        )

        est_stats = _comparison_stats(dip["est"], dip["t2p_est"])
        var_stats = _comparison_stats(dip["var"], dip["t2p_var"])
        assert est_stats["max_abs_diff"] < 1e-5
        assert var_stats["max_abs_diff"] < 1e-5

        py_est_delta = np.abs(dip["est"] - no_dip["est"])
        t2p_est_delta = np.abs(dip["t2p_est"] - no_dip["t2p_est"])
        py_var_delta = np.abs(dip["var"] - no_dip["var"])
        t2p_var_delta = np.abs(dip["t2p_var"] - no_dip["t2p_var"])
        assert t2p_est_delta.mean() > 0.0
        np.testing.assert_allclose(t2p_est_delta, py_est_delta, atol=1e-5)
        np.testing.assert_allclose(t2p_var_delta, py_var_delta, atol=1e-5)


    @pytest.mark.parametrize("aniso1",  [1.0, 0.1])
    @pytest.mark.parametrize("aniso2",  [0.1, 0.01])
    @pytest.mark.parametrize("azimuth", [0.0, 225])
    @pytest.mark.parametrize("dip", [0.0, 5.0])
    def test_ck3d_correctness(self, run_t2p, butte_ck_data, aniso1, aniso2, azimuth, dip):
        """
        3D CK estimates correctness test.
        """

        vgm_p = f"sph 0.0    0.12  {5000*aniso1} 5000 {5000*aniso2} {azimuth} {dip} 0.0 150"   # primary  – Coarse fraction
        vgm_s = f"sph 0.02   0.08  {5000*aniso1} 5000 {5000*aniso2} {azimuth} {dip} 0.0 150"   # secondary – log10(resistivity)
        vgm_x = f"sph 0.0   0.048  {5000*aniso1} 5000 {5000*aniso2} {azimuth} {dip} 0.0"       # cross

        calc_ck_correlation_3d(run_t2p, butte_ck_data, vgm_p=vgm_p, vgm_s=vgm_s, vgm_x=vgm_x)
