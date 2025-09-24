/* MASS_MAPS minimal diagnostics skeleton (clean rewrite) */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>

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

/* Local constants */
#ifndef MASS_MAPS_Z_EPS
#define MASS_MAPS_Z_EPS 1e-8
#endif

/* Forward decls for helpers provided elsewhere */
void set_weight(pos_data *);

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
  /* Replication shift in physical units: integer box counts times InterPartDist */
  double shift[3] = {
      (ii * MyGrids[0].GSglobal[_x_] + subbox.stabl[_x_]) * params.InterPartDist,
      (jj * MyGrids[0].GSglobal[_y_] + subbox.stabl[_y_]) * params.InterPartDist,
      (kk * MyGrids[0].GSglobal[_z_] + subbox.stabl[_z_]) * params.InterPartDist};

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
  /* Precompute chi at segment endpoints and PLC center in physical units */
  double chi_prev_val = ComovingDistance(z_prev);
  double chi_curr_val = ComovingDistance(z_curr);
  double c0_phys = plc.center[0] * params.InterPartDist;
  double c1_phys = plc.center[1] * params.InterPartDist;
  double c2_phys = plc.center[2] * params.InterPartDist;

  /* Determine LPT order compiled in */
  int lpt_order = 1;
#ifdef TWO_LPT
  lpt_order = 2;
#ifdef THREE_LPT
  lpt_order = 3;
#endif
#endif

  /* We'll count local crossings per replication for diagnostics */
  unsigned long long total_crossings_all_reps = 0ULL;

  /* Precompute scale factor F for set_point (1+z at current segment) */
  PRODFLOAT Fseg = (PRODFLOAT)(z_segment + 1.0);

  /* Products-tile traversal is now the only path */
  /*
   * PRODUCTS-TILE PATH: cover all grid elements in this task's FFT tile.
   * Iterate over global lattice coords owned by this rank and read from products[].
   */
  for (int ir = 0; ir < MassMapSegmentReplicationCount; ++ir)
  {
    int rep_id = MassMapSegmentReplications[ir];
    unsigned long long rep_crossings = 0ULL;
    /* Count entries using the user's condition only */

    /* BoxSize per-axis in physical units */
    double Bx = MyGrids[0].GSglobal[_x_] * params.InterPartDist;
    double By = MyGrids[0].GSglobal[_y_] * params.InterPartDist;
    double Bz = MyGrids[0].GSglobal[_z_] * params.InterPartDist;

    int ii = plc.repls[rep_id].i;
    int jj = plc.repls[rep_id].j;
    int kk = plc.repls[rep_id].k;
    double shift[3] = {ii * Bx, jj * By, kk * Bz};

    for (int ig = (int)MyGrids[0].GSstart[_x_]; ig < (int)(MyGrids[0].GSstart[_x_] + MyGrids[0].GSlocal[_x_]); ++ig)
    {
      for (int jg = (int)MyGrids[0].GSstart[_y_]; jg < (int)(MyGrids[0].GSstart[_y_] + MyGrids[0].GSlocal[_y_]); ++jg)
      {
        for (int kg = (int)MyGrids[0].GSstart[_z_]; kg < (int)(MyGrids[0].GSstart[_z_] + MyGrids[0].GSlocal[_z_]); ++kg)
        {
          pos_data obj;
          mass_maps_set_point_from_products_global(ig, jg, kg, Fseg, &obj);
          if (obj.M == 0)
          {
            continue; /* safety */
          }

          double pos_prev[3], pos_curr[3];
          mass_maps_compute_prev_curr_positions(&obj, pos_prev, pos_curr, lpt_order);

          /* Apply replication shift (global positions -> replicated positions) */
          double Pprev[3] = {pos_prev[0] + shift[0], pos_prev[1] + shift[1], pos_prev[2] + shift[2]};
          double Pcurr[3] = {pos_curr[0] + shift[0], pos_curr[1] + shift[1], pos_curr[2] + shift[2]};

          /* Compute r_prev/r_curr in physical units (wrt PLC center) */
          double dx0 = Pprev[0] - c0_phys, dy0 = Pprev[1] - c1_phys, dz0 = Pprev[2] - c2_phys;
          double dx1 = Pcurr[0] - c0_phys, dy1 = Pcurr[1] - c1_phys, dz1 = Pcurr[2] - c2_phys;
          double r_prev = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
          double r_curr = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
          /* User's shell selection: r_prev < chi_prev and r_curr > chi_curr */
          if (!(r_prev < chi_prev_val && r_curr > chi_curr_val))
            continue;

          ++rep_crossings;
        }
      }
    }

    total_crossings_all_reps += rep_crossings;
    if (!ThisTask && internal.verbose_level >= VDIAG)
      printf("[%s] MASS_MAPS(products): segment %d rep %d local entries %llu\n", fdate(), segment_index, rep_id, rep_crossings);
  }

  if (!ThisTask && internal.verbose_level >= VDIAG)
    printf("[%s] MASS_MAPS(products): segment %d total local entries across %d reps: %llu\n",
           fdate(), segment_index, MassMapSegmentReplicationCount, total_crossings_all_reps);
  return; /* products path handled this segment */
}

#endif /* MASS_MAPS */
