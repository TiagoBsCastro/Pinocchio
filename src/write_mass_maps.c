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

/* Global sheet array */
MassSheet *MassSheets = NULL;
int NMassSheets = 0;
double *MassMapBoundaryZ = NULL;
double *MassMapBoundaryChi = NULL;
double *MassMapBoundaryDA = NULL;

/* Epsilon for duplicate detection (user chose 1e-3) */
#define MASS_MAPS_Z_EPS 1e-3

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

#endif /* MASS_MAPS */
