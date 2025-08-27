/* MASS_MAPS minimal diagnostics skeleton (clean rewrite) */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>

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

/* Global sheet array */
MassSheet *MassSheets = NULL;
int NMassSheets = 0;
double *MassMapBoundaryZ = NULL;
double *MassMapBoundaryChi = NULL;
double *MassMapBoundaryDA = NULL;
/* Per-segment replication candidate list (radial overlap with segment chi span) */
static int *MassMapSegmentReplications = NULL;
static int MassMapSegmentReplicationCount = 0;
static int MassMapSegmentReplicationsCapacity = 0;

/* Epsilon for duplicate detection (user chose 1e-3) */
#define MASS_MAPS_Z_EPS 1e-3

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
    const double eps_match = 1e-6; /* strict: outputs are usually exact; relax if needed */
    double z_first = outputs.z[0];
    double z_last = outputs.z[outputs.n - 1];
    if (fabs(z_first - params.StartingzForPLC) > eps_match)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: First output redshift %g does not match StartingzForPLC=%g (|Δ|=%g > %g).\n",
                z_first, params.StartingzForPLC, fabs(z_first - params.StartingzForPLC), eps_match);
      return 1;
    }
    if (fabs(z_last - params.LastzForPLC) > eps_match)
    {
      if (!ThisTask)
        fprintf(stderr, "MASS_MAPS ERROR: Last output redshift %g does not match LastzForPLC=%g (|Δ|=%g > %g).\n",
                z_last, params.LastzForPLC, fabs(z_last - params.LastzForPLC), eps_match);
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

#define MASS_MAPS_F_EPS 1e-7

static inline double plc_cos_aperture(void)
{
  return cos(params.PLCAperture);
}

static inline void plc_axis(double v[3])
{
  v[0] = plc.zvers[0];
  v[1] = plc.zvers[1];
  v[2] = plc.zvers[2];
}

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
int mass_maps_replication_overlaps_sheet(int rep_id, int sheet_id)
{
  if (rep_id < 0 || rep_id >= plc.Nreplications || sheet_id < 0 || sheet_id >= NMassSheets)
    return 0;
  /* New interpretation: F1,F2 store scale factors (1+z). Convert to z then comoving distance. */
  double F1 = plc.repls[rep_id].F1;
  double F2 = plc.repls[rep_id].F2;
  double z1 = F1 - 1.0;
  double z2 = F2 - 1.0;
  if (z1 < 0)
    z1 = 0; /* clamp */
  if (z2 < 0)
    z2 = 0;
  double chi1 = ComovingDistance(z1);
  double chi2 = ComovingDistance(z2);
  double r_rep_min = (chi1 < chi2) ? chi1 : chi2;
  double r_rep_max = (chi1 > chi2) ? chi1 : chi2;
  double chi_lo = MassSheets[sheet_id].chi_lo;
  double chi_hi = MassSheets[sheet_id].chi_hi;
  /* Overlap if intervals intersect */
  if (r_rep_max < chi_lo || r_rep_min > chi_hi)
    return 0;
  return 1;
}

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

/* Compute F = cos(theta) - cos(aperture) for a given position */
/*
 * mass_maps_compute_F
 * -------------------
 * Compute scalar F = cos(theta) - cos(aperture) for position relative to
 * the PLC observer and axis. Inside-cone test is F >= 0.
 * Returns 1.0 (definitely inside) at the cone apex to avoid division by zero.
 */
static inline double mass_maps_compute_F(const double pos[3])
{
  double axis[3];
  plc_axis(axis);
  double dx[3] = {pos[0] - plc.center[0], pos[1] - plc.center[1], pos[2] - plc.center[2]};
  double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  if (r2 == 0.0)
    return 1.0; /* treat at apex as inside */
  double rinv = 1.0 / sqrt(r2);
  double costh = (dx[0] * axis[0] + dx[1] * axis[1] + dx[2] * axis[2]) * rinv;
  return costh - plc_cos_aperture();
}

/*
 * mass_maps_sign
 * --------------
 * Map floating value F to {-1,0,+1} with tolerance MASS_MAPS_F_EPS.
 */
static inline int mass_maps_sign(double F)
{
  if (fabs(F) < MASS_MAPS_F_EPS)
    return 0;
  return (F > 0.0) ? 1 : -1;
}

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
                                   double *alpha_out,
                                   double entry_pos[3])
{
  if (rep_id < 0 || rep_id >= plc.Nreplications)
    return 0;
  int ii = plc.repls[rep_id].i;
  int jj = plc.repls[rep_id].j;
  int kk = plc.repls[rep_id].k;
  double shift[3] = {ii * MyGrids[0].GSglobal[_x_], jj * MyGrids[0].GSglobal[_y_], kk * MyGrids[0].GSglobal[_z_]};

  double pos_prev[3] = {q[0] + disp_prev[0] + shift[0], q[1] + disp_prev[1] + shift[1], q[2] + disp_prev[2] + shift[2]};
  double pos_curr[3] = {q[0] + disp_curr[0] + shift[0], q[1] + disp_curr[1] + shift[1], q[2] + disp_curr[2] + shift[2]};

  double F_prev = mass_maps_compute_F(pos_prev);
  double F_curr = mass_maps_compute_F(pos_curr);
  int s_prev = mass_maps_sign(F_prev);
  int s_curr = mass_maps_sign(F_curr);
  if (!(s_prev < 0 && s_curr >= 0))
    return 0; /* no crossing of interest */
  double dF = F_curr - F_prev;
  if (fabs(dF) < MASS_MAPS_F_EPS)
    return 0;                  /* ambiguous; skip */
  double alpha = -F_prev / dF; /* fraction from prev->curr to crossing */
  if (alpha < 0.0)
    alpha = 0.0;
  else if (alpha > 1.0)
    alpha = 1.0;
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

/* Check pure geometric inside (angle + radial within PLC range) for a single position */
/*
 * mass_maps_point_inside_lightcone
 * --------------------------------
 * Full geometric test for a single Eulerian position:
 *   - Inside angular aperture (F >= 0)
 *   - Radial distance within global PLC radial span [chi_min, chi_max]
 * Does NOT test sheet membership (which is a separate step) or replication.
 * Returns 1 if inside, else 0.
 */
int mass_maps_point_inside_lightcone(const double pos[3], long *ipix_out)
{
  double F = mass_maps_compute_F(pos);
  if (F < 0.0)
    return 0; /* outside aperture */
  /* radial check relative to overall PLC bounds (Fstart/Fstop store scale factors 1+z) */
  double dx = pos[0] - plc.center[0];
  double dy = pos[1] - plc.center[1];
  double dz = pos[2] - plc.center[2];
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  double chi_min = ComovingDistance(params.LastzForPLC);
  double chi_max = ComovingDistance(params.StartingzForPLC);
  if (r < chi_min || r > chi_max)
    return 0;
  if (ipix_out)
  {
    if (params.MassMapNSIDE > 0)
    {
      /* Convert Cartesian to theta,phi (HEALPix uses theta: [0,pi] from +z, phi: [0,2pi)) */
      double vx = dx / r, vy = dy / r, vz = dz / r;
      /* Direction cosines relative to PLC axis frame: we have plc.zvers as axis; need global spherical relative to +Z.
         Simplest: use global coordinates directly. */
      double theta = acos(fmax(-1.0, fmin(1.0, vz))); /* clamp for numerical safety */
      double phi = atan2(vy, vx);
      if (phi < 0)
        phi += 2.0 * M_PI;
      long ipix; /* CHEALPix: choose RING ordering (widely interoperable) */
      ang2pix_ring(params.MassMapNSIDE, theta, phi, &ipix);
      *ipix_out = ipix;
    }
    else
    {
      *ipix_out = -1;
    }
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
#ifdef MASS_MAPS
  if (NMassSheets <= 0 || params.MassMapNSIDE <= 0)
    return; /* feature disabled or not initialized */
  if (is_first_segment)
    return; /* need previous segment for crossings */
  /* STEP 1: select replication candidates overlapping this segment radial span */
  mass_maps_select_segment_replications(segment_index);
  if (!ThisTask)
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
  /* Later steps will loop only over MassMapSegmentReplications[0..MassMapSegmentReplicationCount-1] */
#endif /* MASS_MAPS */
}

#endif /* MASS_MAPS */
