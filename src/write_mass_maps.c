/*
 * MASS_MAPS module (products-only, per-segment PLC entries → HEALPix maps)
 * -----------------------------------------------------------------------
 * Purpose:
 *   Detect entries of mass elements into a light-cone (PLC) between adjacent
 *   fragmentation outputs and accumulate counts on per-segment HEALPix maps.
 *
 * Data flow (per segment):
 *   - Build candidate replications by radial overlap with the segment chi-span.
 *   - Partition the local tile into Lagrangian blocks; for each block+replication,
 *     build a block AABB in physical units expanded by a fixed buffer
 *     (absolute and/or fraction of block half-diagonal) plus MASS_MAPS_CULL_PAD_PHYS.
 *   - Apply conservative culls before touching particles:
 *       • radial shell test via tight AABB distance bounds to the PLC center
 *       • angular aperture test via a cosine upper bound over the AABB
 *   - For passing pairs, traverse particles; test the user crossing condition
 *     H_prev = chi_prev − r_prev, H_curr = chi_curr − r_curr; on success, linearly
 *     interpolate the entry position, gate by aperture, pixelize, and accumulate.
 *   - MPI-reduce per-rep counters and the map to rank 0; write one FITS per segment.
 *
 * Performance notes:
 *   - Lagrangian block culling prunes most work with cheap block-level radial
 *     and angular tests; replication shifts are applied only to passing blocks.
 *   - OpenMP parallelizes over blocks; optional per-thread HEALPix accumulators
 *     reduce atomic contention on the global map.
 *
 * Environment knobs:
 *   MASS_MAPS_BLOCKS[_X|_Y|_Z]: number of Lagrangian blocks per axis (>=1)
 *   MASS_MAPS_CULL_PAD_PHYS:    padding (phys units) added to Lagrangian AABBs (phys units)
 *   MASS_MAPS_DISP_BUFFER_FIXED_PHYS: if >=0, use this fixed buffer (phys units)
 *   MASS_MAPS_DISP_BUFFER_FRAC:       if >=0, per-block buffer = FRAC * (half-diagonal of block) in phys units
 *
 * Accuracy contract:
 *   All refactors preserve the user-specified crossing criterion and output
 *   accumulation semantics. Any optimization must not change results.
 */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifndef PLC
#error "MASS_MAPS requires PLC"
#endif

#ifndef RECOMPUTE_DISPLACEMENTS
#error "MASS_MAPS requires RECOMPUTE_DISPLACEMENTS"
#endif

#ifndef HAVE_CFITSIO
#error "MASS_MAPS requires HAVE_CFITSIO"
#endif
#include <fitsio.h>
#ifdef MASS_MAPS
#include <chealpix.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

/* Local constants */
#ifndef MASS_MAPS_Z_EPS
#define MASS_MAPS_Z_EPS 1e-8
#endif
/* Crossing epsilon for alpha denominator */
#ifndef MASS_MAPS_F_EPS
#define MASS_MAPS_F_EPS 1e-7
#endif
/* Culling defaults  */
#ifndef MASS_MAPS_BLOCKS
#define MASS_MAPS_BLOCKS 8
#endif
#ifndef MASS_MAPS_CULL_PAD_PHYS
#define MASS_MAPS_CULL_PAD_PHYS 0.01
#endif

/* Forward decls for helpers provided elsewhere */
void set_weight(pos_data *);

/* Global sheet array */
MassSheet *MassSheets = NULL;
int NMassSheets = 0;
double *MassMapBoundaryZ = NULL;
double *MassMapBoundaryChi = NULL;
double *MassMapBoundaryDA = NULL;
/* HEALPix maps: one map per mass sheet (segment). */
static int MassMapNSIDE_current = 0;
static long MassMapNPIX = 0;
static double *MassMapSegmentMaps = NULL;  /* length = NMassSheets * MassMapNPIX */
static double *MassMapReduceBuffer = NULL; /* root-only temporary buffer for reductions */
/* Per-segment replication candidate list (radial overlap with segment chi span) */
static int *MassMapSegmentReplications = NULL;
static int MassMapSegmentReplicationCount = 0;
static int MassMapSegmentReplicationsCapacity = 0;

/* Cached PLC center (phys) and sky basis (e1,e2,e3) for pixelization */
static double MassMapPLC_CenterPhys[3] = {0.0, 0.0, 0.0};
static double MassMapPLC_e1[3] = {1.0, 0.0, 0.0};
static double MassMapPLC_e2[3] = {0.0, 1.0, 0.0};
static double MassMapPLC_e3[3] = {0.0, 0.0, 1.0};
static int MassMapPLC_BasisReady = 0;

static inline void mass_maps_cache_plc_basis_and_center(void)
{
  if (MassMapPLC_BasisReady)
    return;
  MassMapPLC_CenterPhys[0] = plc.center[0] * params.InterPartDist;
  MassMapPLC_CenterPhys[1] = plc.center[1] * params.InterPartDist;
  MassMapPLC_CenterPhys[2] = plc.center[2] * params.InterPartDist;
  /* Use PLC basis computed at initialization (orthonormal and unit) */
  MassMapPLC_e1[0] = plc.xvers[0];
  MassMapPLC_e1[1] = plc.xvers[1];
  MassMapPLC_e1[2] = plc.xvers[2];
  MassMapPLC_e2[0] = plc.yvers[0];
  MassMapPLC_e2[1] = plc.yvers[1];
  MassMapPLC_e2[2] = plc.yvers[2];
  MassMapPLC_e3[0] = plc.zvers[0];
  MassMapPLC_e3[1] = plc.zvers[1];
  MassMapPLC_e3[2] = plc.zvers[2];
  MassMapPLC_BasisReady = 1;
}

/* Return 1 if the absolute position lies within the PLC angular aperture */
static inline int mass_maps_entry_inside_aperture(const double pos[3])
{
  mass_maps_cache_plc_basis_and_center();
  double vx = pos[0] - MassMapPLC_CenterPhys[0];
  double vy = pos[1] - MassMapPLC_CenterPhys[1];
  double vz = pos[2] - MassMapPLC_CenterPhys[2];
  double vnorm = sqrt(vx * vx + vy * vy + vz * vz);
  if (vnorm == 0.0)
    return 0;
  double invv = 1.0 / vnorm;
  double vhat_z = (vx * MassMapPLC_e3[0] + vy * MassMapPLC_e3[1] + vz * MassMapPLC_e3[2]) * invv;
  if (vhat_z > 1.0)
    vhat_z = 1.0;
  else if (vhat_z < -1.0)
    vhat_z = -1.0;
  double angle = acos(vhat_z) * (180.0 / acos(-1.0));
  double A = params.PLCAperture;
  if (A > 90.0)
    A = 90.0;
  return angle <= A;
}
/*
 * mass_maps_init_sheets
 * ---------------------
 * Build the array of MassSheet structures from the global output redshift list.
 * Validates:
 *   - At least 2 output redshifts.
 *   - Strictly descending ordering (z[i] > z[i+1]).
 *   - First / last redshifts match PLC StartingzForPLC / LastzForPLC.
 *   - No adjacent duplicates within MASS_MAPS_Z_EPS.
 * Allocates:
 *   - MassSheets (NMassSheets = outputs.n - 1)
 *   - Boundary arrays (Z, Chi, DA) of length outputs.n.
 * Precomputes per-sheet: chi / da bounds, deltas, inverse delta, chi^3 difference.
 * Returns 0 on success, 1 on validation / allocation failure (after emitting message on task 0).
 */
int mass_maps_init_sheets(void)
{
  /* Preconditions: outputs.z[] filled, descending order, outputs.n >= 2 */
  if (outputs.n < 2)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS ERROR: Need at least 2 output redshifts (found %d) to define mass sheets.\n", outputs.n);
    return 1;
  }

  /* Ensure coverage matches PLC redshift range: first output == StartingzForPLC, last output == LastzForPLC */
  {
    double z_first = outputs.z[0];
    double z_last = outputs.z[outputs.n - 1];
    if (fabs(z_first - params.StartingzForPLC) > MASS_MAPS_Z_EPS)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: First output redshift %g does not match StartingzForPLC=%g (|Δ|=%g > %g).\n",
                z_first, params.StartingzForPLC, fabs(z_first - params.StartingzForPLC), MASS_MAPS_Z_EPS);
      return 1;
    }
    if (fabs(z_last - params.LastzForPLC) > MASS_MAPS_Z_EPS)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Last output redshift %g does not match LastzForPLC=%g (|Δ|=%g > %g).\n",
                z_last, params.LastzForPLC, fabs(z_last - params.LastzForPLC), MASS_MAPS_Z_EPS);
      return 1;
    }
  }

  /* Verify strictly descending and no duplicates within eps */
  for (int i = 0; i < outputs.n - 1; ++i)
  {
    double z_hi = outputs.z[i];
    double z_lo = outputs.z[i + 1];
    if (!(z_hi > z_lo))
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Output redshifts must be in strictly descending order (z[%d]=%g, z[%d]=%g).\n", i, z_hi, i + 1, z_lo);
      return 1;
    }
    if (fabs(z_hi - z_lo) < MASS_MAPS_Z_EPS)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Adjacent output redshifts too close (|%g-%g| < %g).\n", z_hi, z_lo, MASS_MAPS_Z_EPS);
      return 1;
    }
  }
  /*
   * Number of mass sheets equals the number of adjacent redshift pairs
   * in the output list: sheets s = 0..(n-2) are bounded by z[s] and z[s+1].
   * We therefore allocate:
   *  - MassSheets: one entry per sheet (outputs.n - 1)
   *  - Boundary arrays (Z/Chi/DA): one value per boundary, length outputs.n
   *    (there are outputs.n boundaries for outputs.n-1 sheets).
   */
  NMassSheets = outputs.n - 1;
  MassSheets = (MassSheet *)malloc(sizeof(MassSheet) * NMassSheets);
  MassMapBoundaryZ = (double *)malloc(sizeof(double) * outputs.n);
  MassMapBoundaryChi = (double *)malloc(sizeof(double) * outputs.n);
  MassMapBoundaryDA = (double *)malloc(sizeof(double) * outputs.n);
  if (!MassSheets)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS ERROR: Cannot allocate MassSheets (%d entries).\n", NMassSheets);
    return 1;
  }
  if (!MassMapBoundaryZ || !MassMapBoundaryChi || !MassMapBoundaryDA)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS ERROR: Cannot allocate boundary arrays (n=%d).\n", outputs.n);
    return 1;
  }

  /* Fill boundary arrays */
  for (int k = 0; k < outputs.n; ++k)
  {
    double z = outputs.z[k];
    MassMapBoundaryZ[k] = z;
    MassMapBoundaryChi[k] = ComovingDistance(z);
    MassMapBoundaryDA[k] = DiameterDistance(z);
  }

  /*
   * For each mass sheet s (bounded by output boundaries s and s+1),
   * cache the upper/lower redshift (z), comoving distance (chi), and
   * angular-diameter distance (da). Compute the spans Δz and Δχ, check
   * that Δχ>0 (monotonic chi with descending z), store its inverse, and
   * precompute χ_hi^3 − χ_lo^3 (useful for volume-weighted quantities).
   */
  for (int s = 0; s < NMassSheets; ++s)
  {
    MassSheet *ms = &MassSheets[s];
    ms->z_hi = MassMapBoundaryZ[s];
    ms->z_lo = MassMapBoundaryZ[s + 1];
    ms->delta_z = ms->z_hi - ms->z_lo;
    ms->chi_hi = MassMapBoundaryChi[s];
    ms->chi_lo = MassMapBoundaryChi[s + 1];
    ms->delta_chi = ms->chi_hi - ms->chi_lo;
    if (!(ms->delta_chi > 0.0))
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Non-positive comoving distance span for sheet %d (chi_hi=%g chi_lo=%g).\n", s, ms->chi_hi, ms->chi_lo);
      return 1;
    }
    ms->inv_dchi = 1.0 / ms->delta_chi;
    ms->da_hi = MassMapBoundaryDA[s];
    ms->da_lo = MassMapBoundaryDA[s + 1];
    double chi_hi3 = ms->chi_hi * ms->chi_hi * ms->chi_hi;
    double chi_lo3 = ms->chi_lo * ms->chi_lo * ms->chi_lo;
    ms->chi3_diff = chi_hi3 - chi_lo3;
  }

  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    printf("[%s] MASS_MAPS: initialized %d mass sheets\n", fdate(), NMassSheets);
    for (int s = 0; s < NMassSheets; ++s)
    {
      MassSheet *ms = &MassSheets[s];
      printf("  sheet %3d: z_hi=%8.4f z_lo=%8.4f Δz=%8.4f chi_hi=%10.4f chi_lo=%10.4f Δchi=%10.4f\n", s, ms->z_hi, ms->z_lo, ms->delta_z, ms->chi_hi, ms->chi_lo, ms->delta_chi);
    }
  }

  return 0;
}

/*
 * mass_maps_free_sheets
 * ---------------------
 * Free MassSheets and the boundary arrays; safe to call multiple times.
 */
void mass_maps_free_sheets(void)
{
  if (MassSheets)
    free(MassSheets);
  MassSheets = NULL;
  NMassSheets = 0;
  if (MassMapBoundaryZ)
    free(MassMapBoundaryZ);
  if (MassMapBoundaryChi)
    free(MassMapBoundaryChi);
  if (MassMapBoundaryDA)
    free(MassMapBoundaryDA);
  MassMapBoundaryZ = MassMapBoundaryChi = MassMapBoundaryDA = NULL;
}

/* ---------------------------------------------------------- */
/* HEALPix map allocation (one map per mass sheet)            */
/* ---------------------------------------------------------- */

/* Return number of pixels for an NSIDE; avoid chealpix dependency in calc */
static inline long mass_maps_npix_from_nside(int nside)
{
  return 12L * (long)nside * (long)nside;
}

/* Accessor: pointer to segment s map (size MassMapNPIX) */
static inline double *mass_maps_segment_ptr(int s)
{
  if (!MassMapSegmentMaps || s < 0 || s >= NMassSheets)
    return NULL;
  return MassMapSegmentMaps + ((long)s) * MassMapNPIX;
}

/* Zero all maps */
static void mass_maps_zero_all_maps(void)
{
  if (!MassMapSegmentMaps)
    return;
  long total = (long)NMassSheets * MassMapNPIX;
  for (long i = 0; i < total; ++i)
    MassMapSegmentMaps[i] = 0.0;
}

/* Free maps */
void mass_maps_free_maps(void)
{
  if (MassMapSegmentMaps)
    free(MassMapSegmentMaps);
  MassMapSegmentMaps = NULL;
  if (MassMapReduceBuffer)
    free(MassMapReduceBuffer);
  MassMapReduceBuffer = NULL;
  MassMapNSIDE_current = 0;
  MassMapNPIX = 0;
}

/* Ensure per-segment HEALPix maps are allocated for given NSIDE */
int mass_maps_init_healpix_maps(int nside)
{
  if (nside <= 0)
    return 1;
  if (NMassSheets <= 0)
    return 1;
  if (MassMapSegmentMaps && MassMapNSIDE_current == nside)
    return 0; /* already ready */

  /* (Re)allocate */
  mass_maps_free_maps();
  MassMapNSIDE_current = nside;
  MassMapNPIX = mass_maps_npix_from_nside(nside);
  long total = (long)NMassSheets * MassMapNPIX;
  MassMapSegmentMaps = (double *)malloc(sizeof(double) * (size_t)total);
  if (!MassMapSegmentMaps)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS ERROR: cannot allocate HEALPix maps (segments=%d npix=%ld total=%ld).\n", NMassSheets, MassMapNPIX, total);
    MassMapNSIDE_current = 0;
    MassMapNPIX = 0;
    return 1;
  }
  mass_maps_zero_all_maps();

  if (!ThisTask && internal.verbose_level >= VDIAG)
    printf("[%s] MASS_MAPS: allocated HEALPix maps NSIDE=%d npix=%ld segments=%d\n", fdate(), nside, MassMapNPIX, NMassSheets);
  return 0;
}

/* ---------------------------------------------------------- */
/* PLC-oriented basis and pixelization helpers                */
/* ---------------------------------------------------------- */

/* Compute HEALPix pixel for segment s at entry between Pprev/Pcurr using alpha in [0,1].
   Orientation is defined by the PLC axis used for groups. Returns 1 on success, 0 on failure. */
/* Debug helpers: compute simple stats and axis pixel id */
static inline void mass_maps_map_stats(const double *map, long npix,
                                       double *sum_out, double *min_out, double *max_out, long *nnz_out)
{
  double s = 0.0;
  double mn = 0.0;
  double mx = 0.0;
  long nnz = 0;
  if (npix > 0)
  {
    mn = mx = map[0];
    for (long i = 0; i < npix; ++i)
    {
      double v = map[i];
      s += v;
      if (v < mn)
        mn = v;
      if (v > mx)
        mx = v;
      if (v != 0.0)
        nnz++;
    }
  }
  if (sum_out)
    *sum_out = s;
  if (min_out)
    *min_out = mn;
  if (max_out)
    *max_out = mx;
  if (nnz_out)
    *nnz_out = nnz;
}

static inline long mass_maps_axis_pixel_ring(int nside)
{
  long pix = -1;
  ang2pix_ring((long)nside, 0.0, 0.0, &pix); /* theta=0 on PLC axis */
  return pix;
}

/* Fast crossing test from absolute, replication-shifted positions with precomputed chi */
static inline int mass_maps_crossing_from_shifted_positions(const double Pprev[3],
                                                            const double Pcurr[3],
                                                            const double c_phys[3],
                                                            double chi_prev,
                                                            double chi_curr,
                                                            double *alpha_out,
                                                            double entry_pos[3])
{
  /* Distances to PLC center */
  double dx0 = Pprev[0] - c_phys[0];
  double dy0 = Pprev[1] - c_phys[1];
  double dz0 = Pprev[2] - c_phys[2];
  double dx1 = Pcurr[0] - c_phys[0];
  double dy1 = Pcurr[1] - c_phys[1];
  double dz1 = Pcurr[2] - c_phys[2];
  double r_prev = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  double r_curr = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  double H_prev = chi_prev - r_prev;
  double H_curr = chi_curr - r_curr;
  if (!(H_prev > 0.0 && H_curr <= 0.0))
    return 0;
  double denom = H_prev - H_curr;
  if (fabs(denom) < MASS_MAPS_F_EPS)
    return 0;
  double alpha = H_prev / denom;
  if (alpha < 0.0)
    alpha = 0.0;
  else if (alpha > 1.0)
    alpha = 1.0;
  if (alpha_out)
    *alpha_out = alpha;
  if (entry_pos)
  {
    entry_pos[0] = Pprev[0] + alpha * (Pcurr[0] - Pprev[0]);
    entry_pos[1] = Pprev[1] + alpha * (Pcurr[1] - Pprev[1]);
    entry_pos[2] = Pprev[2] + alpha * (Pcurr[2] - Pprev[2]);
  }
  return 1;
}

/* ---------------------------------------------------------- */
/* Sub-volume culling helpers (AABB vs PLC shell)             */
/* ---------------------------------------------------------- */

/* Compute tight lower/upper bounds for distance from point c to an AABB [min,max] */
static inline void aabb_distance_bounds_to_point(const double minv[3], const double maxv[3], const double c[3], double *dmin_out, double *dmax_out)
{
  /* d_min: point-to-box distance */
  double d2min = 0.0;
  for (int a = 0; a < 3; ++a)
  {
    double da = 0.0;
    if (c[a] < minv[a])
      da = minv[a] - c[a];
    else if (c[a] > maxv[a])
      da = c[a] - maxv[a];
    d2min += da * da;
  }
  /* d_max: farthest corner distance */
  double d2max = 0.0;
  for (int a = 0; a < 3; ++a)
  {
    double da0 = fabs(minv[a] - c[a]);
    double da1 = fabs(maxv[a] - c[a]);
    double da = (da0 > da1) ? da0 : da1;
    d2max += da * da;
  }
  if (dmin_out)
    *dmin_out = sqrt(d2min);
  if (dmax_out)
    *dmax_out = sqrt(d2max);
}

/* Read integer env var with default fallback; returns >=0 (0 means disabled) */
static inline int getenv_int_default(const char *name, int defval)
{
  const char *v = getenv(name);
  if (!v || !*v)
    return defval;
  int val = atoi(v);
  return val < 0 ? defval : val;
}

static inline double getenv_double_default(const char *name, double defval)
{
  const char *v = getenv(name);
  if (!v || !*v)
    return defval;
  return atof(v);
}

/* Compute HEALPix pixel directly from an absolute position (phys units) */
int mass_maps_compute_pixel_from_pos(const double pos[3], int nside, int nest, long *pix_out)
{
  if (!pix_out || nside <= 0)
    return 0;
  /* Cache center/basis once per process to avoid recomputing in hot paths */
  mass_maps_cache_plc_basis_and_center();
  double vx = pos[0] - MassMapPLC_CenterPhys[0];
  double vy = pos[1] - MassMapPLC_CenterPhys[1];
  double vz = pos[2] - MassMapPLC_CenterPhys[2];
  double vnorm = sqrt(vx * vx + vy * vy + vz * vz);
  if (vnorm == 0.0)
    return 0;
  double invv = 1.0 / vnorm;
  double vhat[3] = {vx * invv, vy * invv, vz * invv};
  const double *e1 = MassMapPLC_e1;
  const double *e2 = MassMapPLC_e2;
  const double *e3 = MassMapPLC_e3;
  double cos_theta = vhat[0] * e3[0] + vhat[1] * e3[1] + vhat[2] * e3[2];
  if (cos_theta > 1.0)
    cos_theta = 1.0;
  else if (cos_theta < -1.0)
    cos_theta = -1.0;
  double theta = acos(cos_theta);
  double x1 = vhat[0] * e1[0] + vhat[1] * e1[1] + vhat[2] * e1[2];
  double y1 = vhat[0] * e2[0] + vhat[1] * e2[1] + vhat[2] * e2[2];
  double phi = atan2(y1, x1);
  if (phi < 0.0)
  {
    double twopi = 2.0 * acos(-1.0);
    phi += twopi;
  }
  long ipix = -1;
  if (nest)
    ang2pix_nest((long)nside, theta, phi, &ipix);
  else
    ang2pix_ring((long)nside, theta, phi, &ipix);
  if (ipix < 0)
    return 0;
  *pix_out = ipix;
  return 1;
}

/* Accumulate a unit weight into the HEALPix map for segment s at the interpolated entry point */
static inline void mass_maps_accumulate_entry_pos(int s,
                                                  const double entry_pos[3],
                                                  double weight)
{
  if (s < 0 || s >= NMassSheets || MassMapNSIDE_current <= 0)
    return;
  long ipix;
  if (!mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, 0 /*RING*/, &ipix))
    return;
  double *map = mass_maps_segment_ptr(s);
  if (!map)
    return;
  map[ipix] += weight;
}

/* Write a single segment's HEALPix map (RING ordering) to a FITS file on task 0 */
int mass_maps_write_segment_map(int s)
{
  if (ThisTask != 0)
    return 0;
  if (s < 0 || s >= NMassSheets)
    return 1;
  if (!MassMapSegmentMaps || MassMapNSIDE_current <= 0 || MassMapNPIX <= 0)
    return 1;

  char fname[3 * LBLENGTH];
  snprintf(fname, sizeof(fname), "pinocchio.%s.massmap.seg%03d.fits", params.RunFlag, s);
  char fitsname[3 * LBLENGTH + 2];
  snprintf(fitsname, sizeof(fitsname), "!%s", fname); /* overwrite if exists */

  fitsfile *fptr = NULL;
  int status = 0;
  long naxes[1];
  naxes[0] = MassMapNPIX;
  double *map = mass_maps_segment_ptr(s);
  if (!map)
    return 1;

  fits_create_file(fptr ? &fptr : &fptr, fitsname, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return 1;
  }
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    return 1;
  }

  /* HEALPix metadata */
  char pixtype[] = "HEALPIX";
  char ordering[] = "RING";
  fits_update_key(fptr, TSTRING, "PIXTYPE", pixtype, "HEALPix pixelisation", &status);
  fits_update_key(fptr, TSTRING, "ORDERING", ordering, "Pixel ordering scheme (RING/NESTED)", &status);
  fits_update_key(fptr, TINT, "NSIDE", &MassMapNSIDE_current, "HEALPix NSIDE", &status);
  long firstpix = 0;
  long lastpix = MassMapNPIX - 1;
  fits_update_key(fptr, TLONG, "FIRSTPIX", &firstpix, "First pixel number (0-based)", &status);
  fits_update_key(fptr, TLONG, "LASTPIX", &lastpix, "Last pixel number (0-based)", &status);

  /* Write data */
  long fpixel = 1; /* FITS 1-based */
  long nelem = MassMapNPIX;
  fits_write_img(fptr, TDOUBLE, fpixel, nelem, map, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    return 1;
  }

  fits_close_file(fptr, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return 1;
  }
  if (internal.verbose_level >= VDIAG)
    printf("[%s] MASS_MAPS: wrote %s (NSIDE=%d NPIX=%ld)\n", fdate(), fname, MassMapNSIDE_current, MassMapNPIX);
  return 0;
}

/*
 * mass_maps_write_sheet_table
 * ---------------------------
 * Rank 0 only: write an ASCII table describing each mass sheet to
 *   pinocchio.<RunFlag>.sheets.out
 * Columns include redshift bounds, comoving / angular diameter distances,
 * and precomputed delta / inverse spans. Intended for reproducibility and
 * external analysis tools.
 * Returns 0 on success (or if nothing to do), 1 on I/O error.
 */
int mass_maps_write_sheet_table(void)
{
  if (ThisTask != 0)
    return 0;
  if (!MassSheets || NMassSheets <= 0)
    return 0;

  char fname[2 * LBLENGTH];
  snprintf(fname, sizeof(fname), "pinocchio.%s.sheets.out", params.RunFlag);
  FILE *fd = fopen(fname, "w");
  if (!fd)
  {
    fprintf(stderr, "Task 0 ERROR: cannot open %s for writing.\n", fname);
    return 1;
  }
  fprintf(fd, "# Mass sheet definitions\n");
  fprintf(fd, "# Columns:\n");
  fprintf(fd, "#  %-3s %-10s %-10s %-10s %-12s %-12s %-12s %-14s %-12s %-12s %-12s\n",
          "id", "z_hi", "z_lo", "delta_z", "chi_hi", "chi_lo", "delta_chi", "inv_delta_chi", "da_hi", "da_lo", "chi3_diff");
  for (int s = 0; s < NMassSheets; ++s)
  {
    MassSheet *ms = &MassSheets[s];
    fprintf(fd, "%3d %-10.8f %-10.8f %-10.8f %-12.8f %-12.8f %-12.8f %-14.12g %-12.8f %-12.8f %-12.8g\n",
            s, ms->z_hi, ms->z_lo, ms->delta_z, ms->chi_hi, ms->chi_lo, ms->delta_chi, ms->inv_dchi, ms->da_hi, ms->da_lo, ms->chi3_diff);
  }
  fclose(fd);
  if (internal.verbose_level >= VDIAG)
    printf("[%s] MASS_MAPS: wrote %s (N=%d)\n", fdate(), fname, NMassSheets);
  return 0;
}

/* ---------------------------------------------------------- */
/* Auxiliary geometry / crossing utilities (minimal version) */
/* ---------------------------------------------------------- */

/* Test if a replication radial interval overlaps sheet radial interval */
/*
 * mass_maps_replication_overlaps_sheet
 * ------------------------------------
 * Quick radial overlap test between a replication volume and a mass sheet.
 * Replication radial interval is reconstructed from plc.repls[rep].F1 / F2
 * (stored with sign convention set during geometry build).
 * Sheet radial interval is [chi_lo, chi_hi].
 * Returns 1 if intervals intersect, else 0. Performs bounds checks.
 * NOTE: This is an approximate radial filter; angular aperture is NOT tested here.
 */
/* (removed) mass_maps_replication_overlaps_sheet: unused after simplifying selection */

/* Return replication radial chi bounds via stored scale factors (1+z) */
static inline void mass_maps_replication_chi_bounds(int rep_id, double *rmin, double *rmax)
{
  double F1 = plc.repls[rep_id].F1;
  double F2 = plc.repls[rep_id].F2;
  double z1 = F1 - 1.0;
  if (z1 < 0)
    z1 = 0;
  double z2 = F2 - 1.0;
  if (z2 < 0)
    z2 = 0;
  double chi1 = ComovingDistance(z1);
  double chi2 = ComovingDistance(z2);
  if (chi1 < chi2)
  {
    *rmin = chi1;
    *rmax = chi2;
  }
  else
  {
    *rmin = chi2;
    *rmax = chi1;
  }
}

/* Build list of replications whose chi interval overlaps this segment chi span */
static void mass_maps_select_segment_replications(int segment_index)
{
  MassMapSegmentReplicationCount = 0;
  if (segment_index <= 0)
    return; /* first segment has no crossings */
  /* Segment spans between previous and current fragmentation redshifts */
  double z_prev = ScaleDep.z[segment_index - 1];
  double z_curr = ScaleDep.z[segment_index];
  /* Convert to chi; ComovingDistance grows with z so chi_prev > chi_curr given descending z array. */
  double chi_prev = ComovingDistance(z_prev);
  double chi_curr = ComovingDistance(z_curr);
  /* Segment radial span is [chi_curr, chi_prev] (chi_prev >= chi_curr). */
  double seg_chi_min = (chi_curr < chi_prev) ? chi_curr : chi_prev;
  double seg_chi_max = (chi_curr > chi_prev) ? chi_curr : chi_prev;
  /* Ensure capacity */
  if (MassMapSegmentReplicationsCapacity < plc.Nreplications)
  {
    int newcap = plc.Nreplications;
    int *tmp = (int *)realloc(MassMapSegmentReplications, sizeof(int) * newcap);
    if (!tmp)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS WARNING: cannot allocate segment replication list (N=%d).\n", plc.Nreplications);
      return;
    }
    MassMapSegmentReplications = tmp;
    MassMapSegmentReplicationsCapacity = newcap;
  }
  for (int r = 0; r < plc.Nreplications; ++r)
  {
    double rmin, rmax;
    mass_maps_replication_chi_bounds(r, &rmin, &rmax);
    if (rmax < seg_chi_min || rmin > seg_chi_max)
      continue; /* no overlap */
    MassMapSegmentReplications[MassMapSegmentReplicationCount++] = r;
  }
  if (!ThisTask && internal.verbose_level >= VDBG)
  {
    printf("[%s] MASS_MAPS: segment %d chi-span=[%g,%g] replication candidates=%d/%d\n",
           fdate(), segment_index, seg_chi_min, seg_chi_max,
           MassMapSegmentReplicationCount, plc.Nreplications);
  }
}

/* ---------------------------------------------------------- */
/* Displacement / position helpers                             */
/* ---------------------------------------------------------- */
/*
 * mass_maps_compute_prev_curr_positions
 * -------------------------------------
 * For a given particle (pos_data), build its Eulerian positions at the
 * previous segment endpoint and at the current segment endpoint, using the
 * stored displacement components (v_prev / v and higher LPT orders) when
 * RECOMPUTE_DISPLACEMENTS is active. These correspond to the two snapshots
 * between which we test for PLC entry. The order argument controls how many
 * LPT orders to include (1,2,3).
 *
 * NOTE: Replication shifts are NOT applied here; caller adds them per
 * replication. Periodic wrapping is also left to the caller if required.
 */
static inline void mass_maps_compute_prev_curr_positions(const pos_data *p,
                                                         double pos_prev[3],
                                                         double pos_curr[3],
                                                         int order)
{
  /* Sum displacements per order (endpoints: weights collapse to selecting
     either *_prev or current components, so no interpolation weights used). */
  for (int i = 0; i < 3; ++i)
  {
    double dp_prev = p->v_prev[i];
    double dp_curr = p->v[i];
#ifdef TWO_LPT
    if (order > 1)
    {
      dp_prev += p->v2_prev[i];
      dp_curr += p->v2[i];
    }
#ifdef THREE_LPT
    if (order > 2)
    {
      dp_prev += p->v31_prev[i] + p->v32_prev[i];
      dp_curr += p->v31[i] + p->v32[i];
    }
#endif /* THREE_LPT */
#endif /* TWO_LPT */
    pos_prev[i] = p->q[i] + dp_prev;
    pos_curr[i] = p->q[i] + dp_curr;
  }
  /* Rescale from lattice units (grid cells) to physical box units using InterPartDist. */
  double scale = params.InterPartDist; /* typical: BoxSize / GridSize */
  pos_prev[0] *= scale;
  pos_prev[1] *= scale;
  pos_prev[2] *= scale;
  pos_curr[0] *= scale;
  pos_curr[1] *= scale;
  pos_curr[2] *= scale;
}

/* ---------------------------------------------------------- */
/* Build pos_data from products (FFT local tile, global coords) */
/* ---------------------------------------------------------- */
/*
 * mass_maps_set_point_from_products_global
 * ----------------------------------------
 * Build a pos_data object for a particle addressed by its GLOBAL lattice
 * coordinates (ig,jg,kg), reading displacement components from products[]
 * at the current task's FFT tile. F is the scale factor (1+z) for the
 * current segment endpoint. q[] is set in GLOBAL lattice units (ig+SHIFT,...).
 * Caller must ensure (ig,jg,kg) lies inside this task's FFT tile.
 */
static inline void mass_maps_set_point_from_products_global(int ig, int jg, int kg,
                                                            PRODFLOAT F,
                                                            pos_data *myobj)
{
  /* Translate global to local FFT indices */
  int li = ig - (int)MyGrids[0].GSstart[_x_];
  int lj = jg - (int)MyGrids[0].GSstart[_y_];
  int lk = kg - (int)MyGrids[0].GSstart[_z_];

  /* Safety: bound check (should not trigger if caller iterates tile properly) */
  if (li < 0 || lj < 0 || lk < 0 ||
      li >= (int)MyGrids[0].GSlocal[_x_] ||
      lj >= (int)MyGrids[0].GSlocal[_y_] ||
      lk >= (int)MyGrids[0].GSlocal[_z_])
  {
    /* Mark as empty */
    myobj->M = 0;
    return;
  }

  unsigned int idx = COORD_TO_INDEX(li, lj, lk, MyGrids[0].GSlocal);

  myobj->z = F - 1.0;
#ifdef SCALE_DEPENDENT
  int S = Smoothing.Nsmooth - 1;
  myobj->myk = Smoothing.k_GM_displ[S];
#else
  myobj->myk = params.k_for_GM;
#endif
  set_weight(myobj);

  myobj->M = 1;
  myobj->q[0] = ig + SHIFT;
  myobj->q[1] = jg + SHIFT;
  myobj->q[2] = kg + SHIFT;
  myobj->v[0] = products[idx].Vel[0];
  myobj->v[1] = products[idx].Vel[1];
  myobj->v[2] = products[idx].Vel[2];
#ifdef TWO_LPT
  myobj->v2[0] = products[idx].Vel_2LPT[0];
  myobj->v2[1] = products[idx].Vel_2LPT[1];
  myobj->v2[2] = products[idx].Vel_2LPT[2];
#ifdef THREE_LPT
  myobj->v31[0] = products[idx].Vel_3LPT_1[0];
  myobj->v31[1] = products[idx].Vel_3LPT_1[1];
  myobj->v31[2] = products[idx].Vel_3LPT_1[2];
  myobj->v32[0] = products[idx].Vel_3LPT_2[0];
  myobj->v32[1] = products[idx].Vel_3LPT_2[1];
  myobj->v32[2] = products[idx].Vel_3LPT_2[2];
#endif
#endif

#ifdef RECOMPUTE_DISPLACEMENTS
  myobj->v_prev[0] = products[idx].Vel_prev[0];
  myobj->v_prev[1] = products[idx].Vel_prev[1];
  myobj->v_prev[2] = products[idx].Vel_prev[2];
#ifdef TWO_LPT
  myobj->v2_prev[0] = products[idx].Vel_2LPT_prev[0];
  myobj->v2_prev[1] = products[idx].Vel_2LPT_prev[1];
  myobj->v2_prev[2] = products[idx].Vel_2LPT_prev[2];
#ifdef THREE_LPT
  myobj->v31_prev[0] = products[idx].Vel_3LPT_1_prev[0];
  myobj->v31_prev[1] = products[idx].Vel_3LPT_1_prev[1];
  myobj->v31_prev[2] = products[idx].Vel_3LPT_1_prev[2];
  myobj->v32_prev[0] = products[idx].Vel_3LPT_2_prev[0];
  myobj->v32_prev[1] = products[idx].Vel_3LPT_2_prev[1];
  myobj->v32_prev[2] = products[idx].Vel_3LPT_2_prev[2];
#endif
#endif
#endif
}

/* ---------------------------------------------------------- */
/* Products-based positions dump (no replication)             */
/* ---------------------------------------------------------- */
void mass_maps_dump_positions_products_for_segment(int s)
{
  if (s <= 0 || s > ScaleDep.myseg)
    return; /* need a previous boundary */

  /* Determine LPT order compiled in */
  int lpt_order = 1;
#ifdef TWO_LPT
  lpt_order = 2;
#ifdef THREE_LPT
  lpt_order = 3;
#endif
#endif

  double z_prev = ScaleDep.z[s - 1];
  double z_curr = ScaleDep.z[s];
  PRODFLOAT Fseg = (PRODFLOAT)(z_curr + 1.0);
  /* For distance diagnostics */
  double chi_prev = ComovingDistance(z_prev);
  double chi_curr = ComovingDistance(z_curr);
  double c0 = plc.center[0] * params.InterPartDist;
  double c1 = plc.center[1] * params.InterPartDist;
  double c2 = plc.center[2] * params.InterPartDist;

  char fname[3 * LBLENGTH];
  snprintf(fname, sizeof(fname), "pinocchio.%s.positions.products.seg%03d.task%05d.out",
           params.RunFlag, s, ThisTask);
  FILE *fd = fopen(fname, "w");
  if (!fd)
  {
    if (!ThisTask)
      fprintf(stderr, "MASS_MAPS dump ERROR: cannot open %s for writing.\n", fname);
    return;
  }
  fprintf(fd, "# products-tile positions dump for segment %d: z_prev=%.8f z_curr=%.8f (Task %d)\n",
          s, z_prev, z_curr, ThisTask);
  /* Box sizes in physical units for diagnostics */
  double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
  double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
  double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;
  fprintf(fd, "# PLC center (phys units): %.9g %.9g %.9g  InterPartDist=%.9g  BoxSize=(%.9g, %.9g, %.9g)\n",
          c0, c1, c2, params.InterPartDist, Bx, By, Bz);
  fprintf(fd, "# ig jg kg  x_prev y_prev z_prev  x_curr y_curr z_curr  r_prev r_curr chi_prev chi_curr\n");

  int ig_start = (int)MyGrids[0].GSstart[_x_];
  int jg_start = (int)MyGrids[0].GSstart[_y_];
  int kg_start = (int)MyGrids[0].GSstart[_z_];
  int ig_end = ig_start + (int)MyGrids[0].GSlocal[_x_];
  int jg_end = jg_start + (int)MyGrids[0].GSlocal[_y_];
  int kg_end = kg_start + (int)MyGrids[0].GSlocal[_z_];

  for (int ig = ig_start; ig < ig_end; ++ig)
  {
    for (int jg = jg_start; jg < jg_end; ++jg)
    {
      for (int kg = kg_start; kg < kg_end; ++kg)
      {
        pos_data obj;
        mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
        if (obj.M == 0)
          continue; /* safety guard */

        double pos_prev[3], pos_curr[3];
        mass_maps_compute_prev_curr_positions(&obj, pos_prev, pos_curr, lpt_order);
        /* Distances to PLC center in physical units (no replication shift here) */
        double dxp0 = pos_prev[0] - c0, dyp0 = pos_prev[1] - c1, dzp0 = pos_prev[2] - c2;
        double dxp1 = pos_curr[0] - c0, dyp1 = pos_curr[1] - c1, dzp1 = pos_curr[2] - c2;
        double r_prev = sqrt(dxp0 * dxp0 + dyp0 * dyp0 + dzp0 * dzp0);
        double r_curr = sqrt(dxp1 * dxp1 + dyp1 * dyp1 + dzp1 * dzp1);

        fprintf(fd, "%d %d %d  %.9g %.9g %.9g  %.9g %.9g %.9g  %.9g %.9g %.9g %.9g\n",
                ig, jg, kg,
                pos_prev[0], pos_prev[1], pos_prev[2],
                pos_curr[0], pos_curr[1], pos_curr[2],
                r_prev, r_curr, chi_prev, chi_curr);
      }
    }
  }
  fclose(fd);
}

/* (removed) angular/radial helper functions: unused with user-based crossing criterion */

/* Determine if particle changed sign (outside->inside) between previous & current displacement states for a replication.
   Inputs are Eulerian displacements already (prev & curr) relative to q (Lagrangian) plus replication shift via rep indices. */
/*
 * mass_maps_particle_sign_change
 * ------------------------------
 * Test whether a particle (given previous & current displacement states) enters
 * the PLC cone during the latest segment for a specific replication.
 *
 * Inputs:
 *   rep_id      - replication index
 *   q[3]        - base (Lagrangian) position of particle
 *   disp_prev   - previous total Eulerian displacement vector
 *   disp_curr   - current total Eulerian displacement vector
 *   alpha_out   - (optional) on success, fraction in [0,1] from previous to current state
 *   entry_pos   - (optional) interpolated Eulerian position at crossing
 *
 * Logic:
 *   1. Build previous and current replicated positions.
 *   2. Compute F_prev, F_curr.
 *   3. Detect sign change (outside -> inside: F_prev < 0, F_curr >= 0).
 *   4. Linear interpolation alpha = -F_prev / (F_curr - F_prev); clamp to [0,1].
 * Skips ambiguous cases where |ΔF| < MASS_MAPS_F_EPS.
 *
 * Returns 1 if a valid crossing detected, else 0.
 */
int mass_maps_particle_sign_change(int rep_id,
                                   const double q[3],
                                   const double disp_prev[3],
                                   const double disp_curr[3],
                                   double z_prev,
                                   double z_curr,
                                   double *alpha_out,
                                   double entry_pos[3])
{
  if (rep_id < 0 || rep_id >= plc.Nreplications)
    return 0;
  int ii = plc.repls[rep_id].i;
  int jj = plc.repls[rep_id].j;
  int kk = plc.repls[rep_id].k;
  /* Replication shift in physical units: integer global box counts only (no subbox offsets) */
  double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
  double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
  double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;
  double shift[3] = {ii * Bx, jj * By, kk * Bz};

  /* If q==0 vector, treat disp_* as absolute positions (already q+dp). Otherwise treat as displacements relative to q. */
  int q_is_zero = (q[0] == 0.0 && q[1] == 0.0 && q[2] == 0.0);
  double pos_prev[3], pos_curr[3];
  if (q_is_zero)
  {
    pos_prev[0] = disp_prev[0] + shift[0];
    pos_prev[1] = disp_prev[1] + shift[1];
    pos_prev[2] = disp_prev[2] + shift[2];
    pos_curr[0] = disp_curr[0] + shift[0];
    pos_curr[1] = disp_curr[1] + shift[1];
    pos_curr[2] = disp_curr[2] + shift[2];
  }
  else
  {
    pos_prev[0] = q[0] + disp_prev[0] + shift[0];
    pos_prev[1] = q[1] + disp_prev[1] + shift[1];
    pos_prev[2] = q[2] + disp_prev[2] + shift[2];
    pos_curr[0] = q[0] + disp_curr[0] + shift[0];
    pos_curr[1] = q[1] + disp_curr[1] + shift[1];
    pos_curr[2] = q[2] + disp_curr[2] + shift[2];
  }

  /* New crossing logic (user condition):
     Define H_prev = chi_prev - r_prev and H_curr = chi_curr - r_curr.
     Crossing occurs if H_prev > 0 and H_curr <= 0. Estimate alpha by
     linear interpolation of H between the two endpoints. */
  double c0 = plc.center[0] * params.InterPartDist;
  double c1 = plc.center[1] * params.InterPartDist;
  double c2 = plc.center[2] * params.InterPartDist;
  double dxp0 = pos_prev[0] - c0, dyp0 = pos_prev[1] - c1, dzp0 = pos_prev[2] - c2;
  double dxp1 = pos_curr[0] - c0, dyp1 = pos_curr[1] - c1, dzp1 = pos_curr[2] - c2;
  double r_prev = sqrt(dxp0 * dxp0 + dyp0 * dyp0 + dzp0 * dzp0);
  double r_curr = sqrt(dxp1 * dxp1 + dyp1 * dyp1 + dzp1 * dzp1);
  double chi_prev = ComovingDistance(z_prev);
  double chi_curr = ComovingDistance(z_curr);
  double H_prev = chi_prev - r_prev;
  double H_curr = chi_curr - r_curr;

  /* Early exit if no crossing according to the user's condition */
  if (!(H_prev > 0.0 && H_curr <= 0.0))
  {
    static int sanity_no_cross_prints = 0;
    if (!ThisTask && internal.verbose_level >= VDIAG && sanity_no_cross_prints < 10)
    {
      printf("MASS_MAPS SANITY(no-entry, user): rep=%d z_prev=%.3f z_curr=%.3f H_prev=% .3e H_curr=% .3e r_prev=%.2f r_curr=%.2f chi_prev=%.2f chi_curr=%.2f\n",
             rep_id, z_prev, z_curr, H_prev, H_curr, r_prev, r_curr, chi_prev, chi_curr);
      sanity_no_cross_prints++;
    }
    return 0;
  }

  double denom = H_prev - H_curr;
  if (fabs(denom) < MASS_MAPS_F_EPS)
    return 0;                    /* ambiguous, skip */
  double alpha = H_prev / denom; /* fraction from prev->curr to crossing */
  /* Clamp to [0,1] just in case */
  if (alpha < 0.0)
    alpha = 0.0;
  else if (alpha > 1.0)
    alpha = 1.0;

  /* Optional sanity print for a few crossings */
  {
    static int sanity_cross_prints = 0;
    if (!ThisTask && internal.verbose_level >= VDIAG && sanity_cross_prints < 10)
    {
      printf("MASS_MAPS SANITY(entry, user): rep=%d z_prev=%.3f z_curr=%.3f H_prev=% .3e H_curr=% .3e alpha=%.3f r_prev=%.2f r_curr=%.2f\n",
             rep_id, z_prev, z_curr, H_prev, H_curr, alpha, r_prev, r_curr);
      sanity_cross_prints++;
    }
  }
  if (alpha_out)
    *alpha_out = alpha;
  if (entry_pos)
  {
    entry_pos[0] = pos_prev[0] + alpha * (pos_curr[0] - pos_prev[0]);
    entry_pos[1] = pos_prev[1] + alpha * (pos_curr[1] - pos_prev[1]);
    entry_pos[2] = pos_prev[2] + alpha * (pos_curr[2] - pos_prev[2]);
  }
  return 1;
}

/* ---------------------------------------------------------- */
/* Orchestrator (skeleton)                                    */
/* ---------------------------------------------------------- */
/*
 * mass_maps_process_segment (SKELETON)
 * -----------------------------------
 * Called once per fragmentation segment immediately AFTER build_groups()
 * when both previous and current displacement states are in memory and
 * group membership for the current segment has been established.
 *
 * Responsibilities (future implementation plan):
 *   1. Lazy allocate / sanity-check mass map accumulation arrays.
 *   2. Build per-particle halo membership mask (in-halo vs field).
 *   3. Determine candidate replication list per sheet (radial prune).
 *   4. Loop particles:
 *        - For each relevant replication, test sign change
 *          (mass_maps_particle_sign_change) for PLC entry.
 *        - If crossing, interpolate entry position, locate sheet
 *          (via MassMapBoundaryChi), compute HEALPix pixel, and
 *          accumulate to total / halo / field counters.
 *   5. (Deferred) Optionally accumulate ancillary stats (e.g. z_entry histogram).
 *   6. Defer MPI reduction & FITS output until all segments processed.
 *
 * Current skeleton: only guards and verbose diagnostic (rank 0).
 */
void mass_maps_process_segment(int segment_index, double z_segment, int is_first_segment)
{
  /* One-time HEALPix map allocator */
  {
    static int hpix_init_done = 0;
    if (!hpix_init_done)
    {
      if (params.MassMapNSIDE > 0)
        mass_maps_init_healpix_maps(params.MassMapNSIDE);
      hpix_init_done = 1;
    }
  }
  /* One-time units diagnostic to stdout */
  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    static int printed_units = 0;
    if (!printed_units)
    {
      double c0 = plc.center[0] * params.InterPartDist;
      double c1 = plc.center[1] * params.InterPartDist;
      double c2 = plc.center[2] * params.InterPartDist;
      double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
      double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
      double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;
      printf("[%s] MASS_MAPS units: PLC center (phys)=(%.9g, %.9g, %.9g)  InterPartDist=%.9g  BoxSize=(%.9g, %.9g, %.9g)\n",
             fdate(), c0, c1, c2, params.InterPartDist, Bx, By, Bz);
      printed_units = 1;
    }
  }

  /* Products-based per-segment positions dump: write once per segment (only at high verbosity) */
  if (internal.verbose_level >= VDIAG)
  {
    static int last_dumped_seg = -1;
    if (segment_index >= 1 && segment_index != last_dumped_seg)
    {
      mass_maps_dump_positions_products_for_segment(segment_index);
      last_dumped_seg = segment_index;
    }
  }

  if (NMassSheets <= 0 || params.MassMapNSIDE <= 0)
    return; /* feature disabled or not initialized */
  if (is_first_segment)
    return; /* need previous segment for crossings */

  /* STEP 1: select replication candidates overlapping this segment radial span */
  mass_maps_select_segment_replications(segment_index);

  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    int show = MassMapSegmentReplicationCount < 10 ? MassMapSegmentReplicationCount : 10;
    printf("[%s] MASS_MAPS: segment %d z=%g replication candidates %d/%d :",
           fdate(), segment_index, z_segment,
           MassMapSegmentReplicationCount, plc.Nreplications);
    for (int i = 0; i < show; ++i)
      printf(" %d", MassMapSegmentReplications[i]);
    if (MassMapSegmentReplicationCount > show)
      printf(" ...");
    printf("\n");
  }

  /* STEP 2: loop over replications and accumulate entries using products over the FFT tile */
  if (MassMapSegmentReplicationCount <= 0)
    return;

  /* Segment endpoint redshifts */
  double z_prev = ScaleDep.z[segment_index - 1];
  double z_curr = z_segment; /* expected to be ScaleDep.z[segment_index] */
  /* Segment chi bounds */
  double chi_prev_seg = ComovingDistance(z_prev);
  double chi_curr_seg = ComovingDistance(z_curr);
  double seg_chi_min = chi_prev_seg < chi_curr_seg ? chi_prev_seg : chi_curr_seg;
  double seg_chi_max = chi_prev_seg > chi_curr_seg ? chi_prev_seg : chi_curr_seg;
  /* PLC center (phys) */
  double c_phys[3] = {plc.center[0] * params.InterPartDist,
                      plc.center[1] * params.InterPartDist,
                      plc.center[2] * params.InterPartDist};

  /* Determine LPT order compiled in */
  int lpt_order = 1;
#ifdef TWO_LPT
  lpt_order = 2;
#ifdef THREE_LPT
  lpt_order = 3;
#endif
#endif

  /* We'll count local crossings per replication and reduce globally */
  unsigned long long total_crossings_all_reps = 0ULL;
  unsigned long long total_crossings_all_reps_global = 0ULL;
  /* Store local per-rep counts to reduce in one shot */
  unsigned long long *rep_crossings_local = NULL;
  unsigned long long *rep_crossings_global = NULL;
  if (MassMapSegmentReplicationCount > 0)
  {
    rep_crossings_local = (unsigned long long *)calloc(MassMapSegmentReplicationCount, sizeof(unsigned long long));
    if (!ThisTask)
      rep_crossings_global = (unsigned long long *)calloc(MassMapSegmentReplicationCount, sizeof(unsigned long long));
  }

  /* Precompute scale factor F for set_point (1+z at current segment) */
  PRODFLOAT Fseg = (PRODFLOAT)(z_segment + 1.0);

  /* Optional: Lagrangian sub-volume culling setup */
  int blocks_all = getenv_int_default("MASS_MAPS_BLOCKS", MASS_MAPS_BLOCKS);
  int blocks_x = getenv_int_default("MASS_MAPS_BLOCKS_X", blocks_all);
  int blocks_y = getenv_int_default("MASS_MAPS_BLOCKS_Y", blocks_all);
  int blocks_z = getenv_int_default("MASS_MAPS_BLOCKS_Z", blocks_all);
  double cull_pad = getenv_double_default("MASS_MAPS_CULL_PAD_PHYS", MASS_MAPS_CULL_PAD_PHYS);
  int thread_accum = getenv_int_default("MASS_MAPS_THREAD_ACCUM", 0); /* per-thread HEALPix accumulators */
  /* Precompute cos(aperture) for conservative angular culling */
  double A = params.PLCAperture;
  if (A > 90.0)
    A = 90.0;
  double cosA = cos(A * (acos(-1.0) / 180.0));
  int nx_local = (int)MyGrids[0].GSlocal[_x_];
  int ny_local = (int)MyGrids[0].GSlocal[_y_];
  int nz_local = (int)MyGrids[0].GSlocal[_z_];
  int gx_start = (int)MyGrids[0].GSstart[_x_];
  int gy_start = (int)MyGrids[0].GSstart[_y_];
  int gz_start = (int)MyGrids[0].GSstart[_z_];
  /* Cap blocks to local sizes */
  if (blocks_x > nx_local)
    blocks_x = nx_local;
  if (blocks_y > ny_local)
    blocks_y = ny_local;
  if (blocks_z > nz_local)
    blocks_z = nz_local;
  int culling_enabled = (blocks_x > 0 && blocks_y > 0 && blocks_z > 0);
  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    printf("[%s] MASS_MAPS culling: blocks=(%d,%d,%d) pad=%.3g enabled=%d (Lagrangian-only)\n",
           fdate(), blocks_x, blocks_y, blocks_z, cull_pad, culling_enabled);
  }

  /* Fixed buffer options: absolute phys value and/or fraction of block size */
  double disp_fixed_phys = getenv_double_default("MASS_MAPS_DISP_BUFFER_FIXED_PHYS", -1.0);
  double disp_frac = getenv_double_default("MASS_MAPS_DISP_BUFFER_FRAC", -1.0);
  /* Last-segment displacement vs fixed-buffer diagnostics */
  double warn_frac = getenv_double_default("MASS_MAPS_DISP_BUFFER_WARN_FRAC", 0.05); /* warn if >= this fraction exceed buffer */
  int last_segment = (segment_index == NMassSheets);
  int collect_disp_stats = (last_segment && disp_fixed_phys > 0.0);
  unsigned long long disp_total_local = 0ULL;  /* number of particles evaluated */
  unsigned long long disp_exceed_local = 0ULL; /* number with |Pcurr-Pprev| > fixed buffer */
  double disp_max_local = 0.0;                 /* maximum |Pcurr-Pprev| observed */
  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    if (disp_fixed_phys >= 0.0)
      printf("[%s] MASS_MAPS culling: using fixed buffer phys=%.6g\n", fdate(), disp_fixed_phys);
    if (disp_frac >= 0.0)
      printf("[%s] MASS_MAPS culling: using fractional buffer frac=%.6g of block half-diagonal\n", fdate(), disp_frac);
    if (collect_disp_stats)
      printf("[%s] MASS_MAPS buffer-diag: last segment will collect displacement stats vs fixed buffer phys=%.6g (warn frac=%.3g)\n",
             fdate(), disp_fixed_phys, warn_frac);
  }
  /* Build block edge indices in local coordinates (inclusive ranges per block) */
  int *edge_x = NULL, *edge_y = NULL, *edge_z = NULL;
  int nbx = 0, nby = 0, nbz = 0;
  if (culling_enabled)
  {
    nbx = blocks_x;
    nby = blocks_y;
    nbz = blocks_z;
    edge_x = (int *)malloc(sizeof(int) * (size_t)(nbx + 1));
    edge_y = (int *)malloc(sizeof(int) * (size_t)(nby + 1));
    edge_z = (int *)malloc(sizeof(int) * (size_t)(nbz + 1));
    if (!edge_x || !edge_y || !edge_z)
    {
      culling_enabled = 0;
    }
    else
    {
      for (int b = 0; b <= nbx; ++b)
        edge_x[b] = (int)((long)b * (long)nx_local / (long)nbx);
      for (int b = 0; b <= nby; ++b)
        edge_y[b] = (int)((long)b * (long)ny_local / (long)nby);
      for (int b = 0; b <= nbz; ++b)
        edge_z[b] = (int)((long)b * (long)nz_local / (long)nbz);
    }
  }

  /* Precompute box shifts for all replication candidates once
   * Rationale: avoids recomputing (ii,jj,kk)*Box for every particle.
   */
  double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
  double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
  double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;
  double *RepShift = NULL;
  if (MassMapSegmentReplicationCount > 0)
  {
    RepShift = (double *)malloc(sizeof(double) * 3 * (size_t)MassMapSegmentReplicationCount);
    if (RepShift)
    {
      for (int ir = 0; ir < MassMapSegmentReplicationCount; ++ir)
      {
        int rep_id = MassMapSegmentReplications[ir];
        int ii = plc.repls[rep_id].i;
        int jj = plc.repls[rep_id].j;
        int kk = plc.repls[rep_id].k;
        RepShift[3 * ir + 0] = ii * Bx;
        RepShift[3 * ir + 1] = jj * By;
        RepShift[3 * ir + 2] = kk * Bz;
      }
    }
  }

  /* Culling-enabled path: iterate blocks first to reuse particle positions across replications
   * Hotspot: the inner loops over particles and replications are dominant.
   * Current approach uses a per-block cache of positions to reduce recomputation
   * and OpenMP over blocks. Map updates default to atomics; an optional per-thread
   * accumulator (MASS_MAPS_THREAD_ACCUM) further reduces contention with a final merge.
   */
  if (culling_enabled && RepShift)
  {
    int s_map = segment_index - 1;
    /* Diagnostics: track block visitation effectiveness */
    long long blocks_total = (long long)nbx * (long long)nby * (long long)nbz;
    long long blocks_init_count = 0;    /* blocks with any particles (AABB initialized) */
    long long blocks_skipped_shell = 0; /* initialized blocks rejected by shell test for all reps */
    long long blocks_visited = 0;       /* initialized blocks with at least one passing replication */
                                        /* Parallelize over blocks; each iteration handles one block fully */
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(dynamic)
#endif
    for (int bz = 0; bz < nbz; ++bz)
      for (int by = 0; by < nby; ++by)
        for (int bx_i = 0; bx_i < nbx; ++bx_i)
        {
          int kg0_local = edge_z[bz];
          int kg1_local = edge_z[bz + 1] - 1;
          if (kg1_local < kg0_local)
            continue;
          int kg0 = gz_start + kg0_local;
          int kg1 = gz_start + kg1_local;
          int jg0_local = edge_y[by];
          int jg1_local = edge_y[by + 1] - 1;
          if (jg1_local < jg0_local)
            continue;
          int jg0 = gy_start + jg0_local;
          int jg1 = gy_start + jg1_local;
          int ig0_local = edge_x[bx_i];
          int ig1_local = edge_x[bx_i + 1] - 1;
          if (ig1_local < ig0_local)
            continue;
          int ig0 = gx_start + ig0_local;
          int ig1 = gx_start + ig1_local;
          /* Determine whether this block has content to consider */
          int block_initialized = 1; /* Lagrangian-only: any non-empty index range is considered initialized */
          if (!block_initialized)
            continue;
#ifdef _OPENMP
#pragma omp atomic update
#endif
          blocks_init_count++;

          /* Build replication pass list for this block */
          int *rep_pass = (int *)alloca(sizeof(int) * (size_t)MassMapSegmentReplicationCount);
          int pass_n = 0;
          for (int ir = 0; ir < MassMapSegmentReplicationCount; ++ir)
          {
            const double *shift = &RepShift[3 * ir];
            /* Lagrangian AABB expanded by fixed buffer and pad (phys units) */
            double scale = params.InterPartDist;
            double min_s[3] = {(ig0 + SHIFT) * scale, (jg0 + SHIFT) * scale, (kg0 + SHIFT) * scale};
            double max_s[3] = {(ig1 + SHIFT) * scale, (jg1 + SHIFT) * scale, (kg1 + SHIFT) * scale};
            /* Compute per-block buffer: fixed phys plus optional fraction of half-diagonal */
            double buf_phys = 0.0;
            if (disp_fixed_phys >= 0.0)
              buf_phys += disp_fixed_phys;
            if (disp_frac >= 0.0)
            {
              double dx = (ig1 - ig0 + 1) * scale * 0.5;
              double dy = (jg1 - jg0 + 1) * scale * 0.5;
              double dz = (kg1 - kg0 + 1) * scale * 0.5;
              double half_diag = sqrt(dx * dx + dy * dy + dz * dz);
              buf_phys += disp_frac * half_diag;
            }
            double pad_total = buf_phys + cull_pad;
            for (int a = 0; a < 3; ++a)
            {
              min_s[a] -= pad_total;
              max_s[a] += pad_total;
              min_s[a] += shift[a];
              max_s[a] += shift[a];
            }
            double dmin, dmax;
            aabb_distance_bounds_to_point(min_s, max_s, c_phys, &dmin, &dmax);
            if (dmin > seg_chi_max || dmax < seg_chi_min)
              continue; /* radial shell rejects */
            /* Conservative angular cull using cosine upper bound:
               If max_{box} (u·r/|r|) ≤ cosA, the whole AABB is outside aperture.
               Upper bound <= (max_u_dot) / (min_norm). */
            double u[3] = {plc.zvers[0], plc.zvers[1], plc.zvers[2]}; /* unit axis */
            /* Build relative bounds r = x - c */
            double rmin[3] = {min_s[0] - c_phys[0], min_s[1] - c_phys[1], min_s[2] - c_phys[2]};
            double rmax[3] = {max_s[0] - c_phys[0], max_s[1] - c_phys[1], max_s[2] - c_phys[2]};
            /* max_u_dot: choose per-axis extreme based on sign of u */
            double max_u_dot = 0.0;
            for (int a = 0; a < 3; ++a)
            {
              double ra = (u[a] >= 0.0) ? rmax[a] : rmin[a];
              max_u_dot += u[a] * ra;
            }
            /* min_norm = dmin from center to AABB; if center inside (dmin==0), skip angular culling */
            int angular_pass = 1;
            if (dmin > 0.0)
            {
              double cos_ub = max_u_dot / dmin;
              if (cos_ub <= cosA)
                angular_pass = 0; /* definitely outside aperture */
            }
            if (!angular_pass)
              continue;
            rep_pass[pass_n++] = ir;
          }
          if (pass_n == 0)
          {
#ifdef _OPENMP
#pragma omp atomic update
#endif
            blocks_skipped_shell++;
            continue; /* no replication needs this block */
          }
#ifdef _OPENMP
#pragma omp atomic update
#endif
          blocks_visited++;

          /* Prepare buffers for this block's particles (thread-local) */
          long nx = (long)(ig1 - ig0 + 1);
          long ny = (long)(jg1 - jg0 + 1);
          long nz = (long)(kg1 - kg0 + 1);
          long nblock = nx * ny * nz;
          if (nblock <= 0)
            continue;
          double *PosPrevBuf = (double *)malloc(sizeof(double) * 3 * (size_t)nblock);
          double *PosCurrBuf = (double *)malloc(sizeof(double) * 3 * (size_t)nblock);
          /* Optional per-thread accumulator for this block to reduce atomics */
          double *LocalMap = NULL;
          if (thread_accum && MassMapNSIDE_current > 0)
            LocalMap = (double *)calloc((size_t)MassMapNPIX, sizeof(double));
          int use_cache = (PosPrevBuf && PosCurrBuf);
          if (use_cache)
          {
            /* Fill cached positions once */
            long idx = 0;
            /* Per-block diagnostics */
            unsigned long long blk_total = 0ULL, blk_exceed = 0ULL;
            double blk_max = 0.0;
            for (int ig = ig0; ig <= ig1; ++ig)
            {
              for (int jg = jg0; jg <= jg1; ++jg)
              {
                for (int kg = kg0; kg <= kg1; ++kg, ++idx)
                {
                  pos_data obj;
                  mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
                  if (obj.M == 0)
                  {
                    /* Mark as zero-length by duplicating prev/curr */
                    PosPrevBuf[3 * idx + 0] = PosPrevBuf[3 * idx + 1] = PosPrevBuf[3 * idx + 2] = 0.0;
                    PosCurrBuf[3 * idx + 0] = PosCurrBuf[3 * idx + 1] = PosCurrBuf[3 * idx + 2] = 0.0;
                    continue;
                  }
                  double pprev[3], pcurr[3];
                  mass_maps_compute_prev_curr_positions(&obj, pprev, pcurr, lpt_order);
                  PosPrevBuf[3 * idx + 0] = pprev[0];
                  PosPrevBuf[3 * idx + 1] = pprev[1];
                  PosPrevBuf[3 * idx + 2] = pprev[2];
                  PosCurrBuf[3 * idx + 0] = pcurr[0];
                  PosCurrBuf[3 * idx + 1] = pcurr[1];
                  PosCurrBuf[3 * idx + 2] = pcurr[2];
                  if (collect_disp_stats)
                  {
                    double dx = pcurr[0] - pprev[0];
                    double dy = pcurr[1] - pprev[1];
                    double dz = pcurr[2] - pprev[2];
                    double d = sqrt(dx * dx + dy * dy + dz * dz);
                    blk_total++;
                    if (d > disp_fixed_phys)
                      blk_exceed++;
                    if (d > blk_max)
                      blk_max = d;
                  }
                }
              }
            }
            if (collect_disp_stats)
            {
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_total_local += blk_total;
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_exceed_local += blk_exceed;
#ifdef _OPENMP
#pragma omp critical
#endif
              {
                if (blk_max > disp_max_local)
                  disp_max_local = blk_max;
              }
            }

            /* Use cached positions for each passing replication */
            for (int ip = 0; ip < pass_n; ++ip)
            {
              int ir = rep_pass[ip];
              const double *shift = &RepShift[3 * ir];
              unsigned long long local_cross = 0ULL;
              long idx2 = 0;
              for (int ig = ig0; ig <= ig1; ++ig)
              {
                for (int jg = jg0; jg <= jg1; ++jg)
                {
                  for (int kg = kg0; kg <= kg1; ++kg, ++idx2)
                  {
                    double Pprev[3] = {PosPrevBuf[3 * idx2 + 0] + shift[0],
                                       PosPrevBuf[3 * idx2 + 1] + shift[1],
                                       PosPrevBuf[3 * idx2 + 2] + shift[2]};
                    double Pcurr[3] = {PosCurrBuf[3 * idx2 + 0] + shift[0],
                                       PosCurrBuf[3 * idx2 + 1] + shift[1],
                                       PosCurrBuf[3 * idx2 + 2] + shift[2]};
                    /* Skip if empty marked */
                    if (Pprev[0] == 0.0 && Pprev[1] == 0.0 && Pprev[2] == 0.0 &&
                        Pcurr[0] == 0.0 && Pcurr[1] == 0.0 && Pcurr[2] == 0.0)
                      continue;
                    double entry_pos[3];
                    if (!mass_maps_crossing_from_shifted_positions(Pprev, Pcurr, c_phys, chi_prev_seg, chi_curr_seg, NULL, entry_pos))
                      continue;
                    /* Angular aperture: only accumulate entries inside PLC cone */
                    if (!mass_maps_entry_inside_aperture(entry_pos))
                      continue;
                    local_cross++;
                    if (s_map >= 0 && s_map < NMassSheets && MassMapNSIDE_current > 0)
                    {
                      long ipix;
                      if (mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, 0, &ipix))
                      {
                        if (LocalMap)
                        {
                          LocalMap[ipix] += 1.0;
                        }
                        else
                        {
                          double *map = mass_maps_segment_ptr(s_map);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                          map[ipix] += 1.0;
                        }
                      }
                    }
                  }
                }
              }
          /* Thread-safe global counters */
#ifdef _OPENMP
#pragma omp atomic update
#endif
              total_crossings_all_reps += local_cross;
              if (rep_crossings_local)
              {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                rep_crossings_local[ir] += local_cross;
              }
            }
          }
          else
          {
            /* No cache available: compute per particle per passing replication */
            /* Per-block diagnostics */
            unsigned long long blk_total = 0ULL, blk_exceed = 0ULL;
            double blk_max = 0.0;
            for (int ip = 0; ip < pass_n; ++ip)
            {
              int ir = rep_pass[ip];
              const double *shift = &RepShift[3 * ir];
              unsigned long long local_cross = 0ULL;
              for (int ig = ig0; ig <= ig1; ++ig)
              {
                for (int jg = jg0; jg <= jg1; ++jg)
                {
                  for (int kg = kg0; kg <= kg1; ++kg)
                  {
                    pos_data obj;
                    mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
                    if (obj.M == 0)
                      continue;
                    double pos_prev[3], pos_curr[3];
                    mass_maps_compute_prev_curr_positions(&obj, pos_prev, pos_curr, lpt_order);
                    if (collect_disp_stats && ip == 0)
                    {
                      double dx = pos_curr[0] - pos_prev[0];
                      double dy = pos_curr[1] - pos_prev[1];
                      double dz = pos_curr[2] - pos_prev[2];
                      double d = sqrt(dx * dx + dy * dy + dz * dz);
                      blk_total++;
                      if (d > disp_fixed_phys)
                        blk_exceed++;
                      if (d > blk_max)
                        blk_max = d;
                    }
                    double Pprev[3] = {pos_prev[0] + shift[0], pos_prev[1] + shift[1], pos_prev[2] + shift[2]};
                    double Pcurr[3] = {pos_curr[0] + shift[0], pos_curr[1] + shift[1], pos_curr[2] + shift[2]};
                    double entry_pos[3];
                    if (!mass_maps_crossing_from_shifted_positions(Pprev, Pcurr, c_phys, chi_prev_seg, chi_curr_seg, NULL, entry_pos))
                      continue;
                    /* Angular aperture: only accumulate entries inside PLC cone */
                    if (!mass_maps_entry_inside_aperture(entry_pos))
                      continue;
                    local_cross++;
                    if (s_map >= 0 && s_map < NMassSheets && MassMapNSIDE_current > 0)
                    {
                      long ipix;
                      if (mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, 0, &ipix))
                      {
                        if (LocalMap)
                        {
                          LocalMap[ipix] += 1.0;
                        }
                        else
                        {
                          double *map = mass_maps_segment_ptr(s_map);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                          map[ipix] += 1.0;
                        }
                      }
                    }
                  }
                }
              }
          /* Thread-safe global counters */
#ifdef _OPENMP
#pragma omp atomic update
#endif
              total_crossings_all_reps += local_cross;
              if (rep_crossings_local)
              {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                rep_crossings_local[ir] += local_cross;
              }
            }
            if (collect_disp_stats)
            {
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_total_local += blk_total;
#ifdef _OPENMP
#pragma omp atomic update
#endif
              disp_exceed_local += blk_exceed;
#ifdef _OPENMP
#pragma omp critical
#endif
              {
                if (blk_max > disp_max_local)
                  disp_max_local = blk_max;
              }
            }
          }
          /* If we used a local map, merge it into the global segment map once per block */
          if (LocalMap)
          {
            double *map = mass_maps_segment_ptr(s_map);
            if (map)
            {
              for (long i = 0; i < MassMapNPIX; ++i)
              {
                double v = LocalMap[i];
                if (v == 0.0)
                  continue;
#ifdef _OPENMP
#pragma omp atomic update
#endif
                map[i] += v;
              }
            }
          }
          if (PosPrevBuf)
            free(PosPrevBuf);
          if (PosCurrBuf)
            free(PosCurrBuf);
          if (LocalMap)
            free(LocalMap);
        }

    /* Print diagnostics on block visitation */
    if (!ThisTask && internal.verbose_level >= VDIAG)
    {
      double init_pct = (blocks_total > 0) ? (100.0 * (double)blocks_init_count / (double)blocks_total) : 0.0;
      double visit_pct = (blocks_total > 0) ? (100.0 * (double)blocks_visited / (double)blocks_total) : 0.0;
      double skip_pct = (blocks_total > 0) ? (100.0 * (double)blocks_skipped_shell / (double)blocks_total) : 0.0;
      printf("[%s] MASS_MAPS culling: blocks total=%lld init=%lld (%.1f%%) visited=%lld (%.1f%%) skipped_by_shell=%lld (%.1f%%)\n",
             fdate(), blocks_total, blocks_init_count, init_pct, blocks_visited, visit_pct, blocks_skipped_shell, skip_pct);
    }
  }
  else
  {
    /* Previous per-rep traversal (fallback without culling or on allocation failure) */
    for (int ir = 0; ir < MassMapSegmentReplicationCount; ++ir)
    {
      int rep_id = MassMapSegmentReplications[ir];
      unsigned long long rep_crossings = 0ULL;
      double shift[3] = {RepShift ? RepShift[3 * ir + 0] : 0.0,
                         RepShift ? RepShift[3 * ir + 1] : 0.0,
                         RepShift ? RepShift[3 * ir + 2] : 0.0};
      unsigned long long local_total = 0ULL, local_exceed = 0ULL;
      double local_max = 0.0;
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static) reduction(+ : rep_crossings, local_total, local_exceed) reduction(max : local_max)
#endif
      for (int ig = gx_start; ig < gx_start + nx_local; ++ig)
      {
        for (int jg = gy_start; jg < gy_start + ny_local; ++jg)
        {
          for (int kg = gz_start; kg < gz_start + nz_local; ++kg)
          {
            pos_data obj;
            mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
            if (obj.M == 0)
              continue; /* safety */
            double pos_prev[3], pos_curr[3];
            mass_maps_compute_prev_curr_positions(&obj, pos_prev, pos_curr, lpt_order);
            if (collect_disp_stats && ir == 0)
            {
              double dx = pos_curr[0] - pos_prev[0];
              double dy = pos_curr[1] - pos_prev[1];
              double dz = pos_curr[2] - pos_prev[2];
              double d = sqrt(dx * dx + dy * dy + dz * dz);
              local_total++;
              if (d > disp_fixed_phys)
                local_exceed++;
              if (d > local_max)
                local_max = d;
            }
            double Pprev[3] = {pos_prev[0] + shift[0], pos_prev[1] + shift[1], pos_prev[2] + shift[2]};
            double Pcurr[3] = {pos_curr[0] + shift[0], pos_curr[1] + shift[1], pos_curr[2] + shift[2]};
            double entry_pos[3];
            if (!mass_maps_crossing_from_shifted_positions(Pprev, Pcurr, c_phys, chi_prev_seg, chi_curr_seg, NULL, entry_pos))
              continue;
            if (!mass_maps_entry_inside_aperture(entry_pos))
              continue;
            ++rep_crossings;
            int s_map = segment_index - 1;
            if (s_map >= 0 && s_map < NMassSheets && MassMapNSIDE_current > 0)
            {
              long ipix;
              if (mass_maps_compute_pixel_from_pos(entry_pos, MassMapNSIDE_current, 0, &ipix))
              {
                double *map = mass_maps_segment_ptr(s_map);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                map[ipix] += 1.0;
              }
            }
          }
        }
      }
      total_crossings_all_reps += rep_crossings;
      if (rep_crossings_local)
        rep_crossings_local[ir] = rep_crossings;
      if (collect_disp_stats && ir == 0)
      {
        disp_total_local += local_total;
        disp_exceed_local += local_exceed;
        if (local_max > disp_max_local)
          disp_max_local = local_max;
      }
      if (!ThisTask && internal.verbose_level >= VDIAG)
      {
        printf("[%s] MASS_MAPS(products): segment %d rep %d local entries %llu\n", fdate(), segment_index, rep_id, rep_crossings);
      }
    }
  }

  if (RepShift)
    free(RepShift);
  /* Free culling structures */
  if (edge_x)
    free(edge_x);
  if (edge_y)
    free(edge_y);
  if (edge_z)
    free(edge_z);

  /* Reduce per-rep and total counts across tasks */
  if (MassMapSegmentReplicationCount > 0)
  {
    /* Per-rep reductions (only if allocations succeeded) */
    if (rep_crossings_local)
    {
      MPI_Reduce(rep_crossings_local, rep_crossings_global,
                 MassMapSegmentReplicationCount, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    /* Total across reps */
    unsigned long long local_total = total_crossings_all_reps;
    MPI_Reduce(&local_total, &total_crossings_all_reps_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (!ThisTask && internal.verbose_level >= VDIAG)
  {
    printf("[%s] MASS_MAPS(products): segment %d total local entries across %d reps: %llu\n",
           fdate(), segment_index, MassMapSegmentReplicationCount, total_crossings_all_reps);
    printf("[%s] MASS_MAPS(products): segment %d total global entries across %d reps: %llu\n",
           fdate(), segment_index, MassMapSegmentReplicationCount, total_crossings_all_reps_global);
  }

  /* Report last-segment displacement vs fixed-buffer stats */
  if (collect_disp_stats)
  {
    unsigned long long disp_total_global = 0ULL, disp_exceed_global = 0ULL;
    double disp_max_global = 0.0;
    MPI_Reduce(&disp_total_local, &disp_total_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disp_exceed_local, &disp_exceed_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disp_max_local, &disp_max_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!ThisTask)
    {
      double frac = (disp_total_global > 0ULL) ? ((double)disp_exceed_global / (double)disp_total_global) : 0.0;
      printf("[%s] MASS_MAPS buffer-diag: last segment displacement vs fixed buffer=%.6g: total=%llu exceed=%llu (%.3f%%) max=%.6g\n",
             fdate(), disp_fixed_phys, disp_total_global, disp_exceed_global, 100.0 * frac, disp_max_global);
      if (disp_total_global > 0ULL && frac >= warn_frac)
      {
        fprintf(stderr, "[%s] MASS_MAPS WARNING: %.2f%% of particles exceeded MASS_MAPS_DISP_BUFFER_FIXED_PHYS=%.6g. Consider increasing the buffer (warn threshold=%.2f%%).\n",
                fdate(), 100.0 * frac, disp_fixed_phys, 100.0 * warn_frac);
      }
    }
  }

  /* Reduce the segment map to root (task 0) */
  if (MassMapNSIDE_current > 0)
  {
    int s_map = segment_index - 1;
    if (s_map >= 0 && s_map < NMassSheets)
    {
      double *map_local = mass_maps_segment_ptr(s_map);
      if (!ThisTask)
      {
        if (!MassMapReduceBuffer)
          MassMapReduceBuffer = (double *)malloc(sizeof(double) * (size_t)MassMapNPIX);
      }
      MPI_Reduce(map_local,
                 (!ThisTask ? MassMapReduceBuffer : NULL),
                 (int)MassMapNPIX,
                 MPI_DOUBLE,
                 MPI_SUM,
                 0,
                 MPI_COMM_WORLD);
      if (!ThisTask && MassMapReduceBuffer)
      {
        /* Copy reduced result back into the segment map on root */
        memcpy(map_local, MassMapReduceBuffer, sizeof(double) * (size_t)MassMapNPIX);
        /* Write to FITS */
        mass_maps_write_segment_map(s_map);
        if (internal.verbose_level >= VDIAG)
        {
          double sum, vmin, vmax;
          long nnz;
          mass_maps_map_stats(map_local, MassMapNPIX, &sum, &vmin, &vmax, &nnz);
          long axis_pix = mass_maps_axis_pixel_ring(MassMapNSIDE_current);
          double axis_val = (axis_pix >= 0 && axis_pix < MassMapNPIX) ? map_local[axis_pix] : -1.0;
          printf("[%s] MASS_MAPS(diag): seg %d map stats: sum=%.0f nnz=%ld min=%.0f max=%.0f axis_pix=%ld axis_val=%.0f\n",
                 fdate(), s_map, sum, nnz, vmin, vmax, axis_pix, axis_val);
        }
      }
    }
  }

  if (rep_crossings_local)
    free(rep_crossings_local);
  if (rep_crossings_global)
    free(rep_crossings_global);
  return; /* products path handled this segment */
}

#endif /* MASS_MAPS */
