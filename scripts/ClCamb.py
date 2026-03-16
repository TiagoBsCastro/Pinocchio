#!/usr/bin/env python3
import argparse
import numpy as np
import camb
from camb import model, sources

C_LIGHT = 299792.458  # km/s

def make_window(z, window_type, results):
    """
    Return normalized W(z) on the provided z-grid, for:
      - 'z'      : flat in redshift (top-hat in z)
      - 'chi'    : flat in comoving distance (W ∝ dχ/dz)
      - 'volume' : flat in comoving volume per sr (W ∝ χ^2 dχ/dz)
    Normalization is ∫dz W(z) = 1, as expected by CAMB number-counts windows.
    """
    # Background quantities
    Hz = results.hubble_parameter(z)     # km/s/Mpc
    chi = results.comoving_radial_distance(z)  # Mpc

    if window_type == 'z':
        W = np.ones_like(z)
    elif window_type == 'chi':
        W = C_LIGHT / Hz                 # ∝ dχ/dz
    elif window_type == 'volume':
        W = chi**2 * (C_LIGHT / Hz)      # ∝ χ^2 dχ/dz
    else:
        raise ValueError("window_type must be one of: 'z', 'chi', 'volume'.")

    # Normalize to unit integral over z
    W = np.array(W, dtype=float)
    W /= np.trapezoid(W, z)
    return W

def build_cls_shell(H0=67.66, Ob0=0.0489, Ocdm0=0.2621, ns=0.9649, As=2.1e-9,
                    mnu=0.06, Neff=3.046, tau=0.054,
                    zmin=0.2, zmax=0.3, nz=256,
                    lmax=2000, limber=True, nonlinear='both',
                    raw_cl=True, window_type='chi'):
    """
    Compute C_ell of the projected matter overdensity shell using a non-flat window.
    window_type ∈ {'z','chi','volume'} as defined above.
    """

    # Convert Ω's to physical densities
    h = H0 / 100.0
    ombh2  = Ob0   * h**2
    omch2  = Ocdm0 * h**2

    # Cosmology and power
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0,
                       ombh2=ombh2,
                       omch2=omch2,
                       mnu=mnu,
                       nnu=Neff,
                       num_massive_neutrinos=1,
                       tau=tau,
                       omk=0.0)
    pars.InitPower.set_params(As=As, ns=ns)

    pars.max_l = int(lmax)
    pars.SourceTerms.limber_windows = bool(limber)
    pars.SourceTerms.limber_phi_lmin = 20

    if nonlinear == 'none':
        pars.NonLinear = model.NonLinear_none
    else:
        pars.NonLinear = model.NonLinear_both

    # Keep only the density contribution; bias=1 -> matter
    st = pars.SourceTerms
    st.counts_density   = True
    st.counts_redshift  = False
    st.counts_lensing   = False
    st.counts_velocity  = False
    st.counts_radial    = False
    st.counts_timedelay = False
    st.counts_ISW       = False
    st.counts_potential = False
    st.counts_evolve    = False

    # Build z-grid and window
    z = np.linspace(zmin, zmax, int(nz))
    # First pass results to get background H(z) and χ(z)
    bg_results = camb.get_results(pars)
    W = make_window(z, window_type, bg_results)

    # Define the counts window with bias=1
    win = sources.SplinedSourceWindow(source_type='counts', bias=1.0)
    win.set_table(z, W)

    # Attach and compute source Cls
    pars.SourceWindows = [win]
    pars.Want_cl_2D_array = True
    results = camb.get_results(pars)  # second call includes sources
    cls_dict = results.get_source_cls_dict(lmax=pars.max_l, raw_cl=bool(raw_cl))

    Cl = cls_dict['W1xW1']
    ells = np.arange(Cl.size)
    return ells, Cl

def main():
    p = argparse.ArgumentParser(description="CAMB: matter shell C_ell with non-flat radial window.")
    # Cosmology defaults you provided
    p.add_argument("--H0",    type=float, default=67.66)
    p.add_argument("--Ob0",   type=float, default=0.0489)
    p.add_argument("--Ocdm0", type=float, default=0.2621)
    p.add_argument("--ns",    type=float, default=0.9649)
    p.add_argument("--As",    type=float, default=2.1e-9)
    p.add_argument("--mnu",   type=float, default=0.06)
    p.add_argument("--Neff",  type=float, default=3.046)
    p.add_argument("--tau",   type=float, default=0.054)

    # Shell and numerics
    p.add_argument("--zmin",  type=float, default=0.2)
    p.add_argument("--zmax",  type=float, default=0.3)
    p.add_argument("--nz",    type=int,   default=256)
    p.add_argument("--lmax",  type=int,   default=2000)
    p.add_argument("--limber", action="store_true")
    p.add_argument("--nonlinear", choices=["none","both"], default="both")
    p.add_argument("--raw_cl", action="store_true")

    # Window choice: z, chi, or volume
    p.add_argument("--window", choices=["z","chi","volume"], default="chi",
                   help="Radial weighting: flat in z, flat in comoving distance, or flat in comoving volume.")

    p.add_argument("--save", type=str, default="")
    args = p.parse_args()

    ells, Cl = build_cls_shell(H0=args.H0, Ob0=args.Ob0, Ocdm0=args.Ocdm0,
                               ns=args.ns, As=args.As, mnu=args.mnu, Neff=args.Neff, tau=args.tau,
                               zmin=args.zmin, zmax=args.zmax, nz=args.nz,
                               lmax=args.lmax, limber=args.limber,
                               nonlinear=args.nonlinear, raw_cl=args.raw_cl,
                               window_type=args.window)

    if args.save:
        np.savetxt(args.save, np.column_stack([ells, Cl]), header="ell  Cl")
        print(f"Saved: {args.save}")
    else:
        print("ell   Cl")
        for L, c in list(zip(ells, Cl))[:10]:
            print(f"{L:4d}  {c:.6e}")

if __name__ == "__main__":
    main()

