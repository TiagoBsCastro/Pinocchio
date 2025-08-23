/* MASS_MAPS minimal diagnostics skeleton (clean rewrite) */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>

#ifndef PLC
#error "MASS_MAPS requires PLC"
#endif

/* ------------------------------------------------------------------------- */
/* Full lattice (products[]) output (production).
  We always attempt to write ALL particles from the initial lattice using
  the products[] ordering because PLC mass map reconstruction requires the
  complete set (ghost layers included). If products[] is unavailable at the
  time of the call (unexpected late free), we fall back to the owned-only
  fragmentation view with a warning. */

/* ------------------------- FULL LATTICE (products[]) --------------------- */
static inline int mass_maps_full_available(void)
{
  return (products != NULL && MyGrids[0].total_local_size > 0);
}

static inline void mass_maps_compute_from_products(int idx,
                                                   unsigned long long *pid,
                                                   double q[3], double x[3])
{
  int ibox, jbox, kbox;
  INDEX_TO_COORD(idx, ibox, jbox, kbox, MyGrids[0].GSlocal);
  const int Gx = MyGrids[0].GSglobal[_x_];
  const int Gy = MyGrids[0].GSglobal[_y_];
  const int Gz = MyGrids[0].GSglobal[_z_];
  const double dx = params.InterPartDist;
  int gx = ibox + MyGrids[0].GSstart[_x_];
  int gy = jbox + MyGrids[0].GSstart[_y_];
  int gz = kbox + MyGrids[0].GSstart[_z_];
  *pid = ((unsigned long long)gx * Gy + gy) * Gz + gz + 1ULL;
  q[0] = gx * dx;
  q[1] = gy * dx;
  q[2] = gz * dx;
  for (int d = 0; d < 3; d++)
  {
    int g = (d == 0 ? gx : (d == 1 ? gy : gz));
    double pos = g + SHIFT + products[idx].Vel[d];
#ifdef TWO_LPT
    pos += products[idx].Vel_2LPT[d];
#ifdef THREE_LPT
    pos += products[idx].Vel_3LPT_1[d] + products[idx].Vel_3LPT_2[d];
#endif
#endif
    int GL = (d == 0 ? Gx : (d == 1 ? Gy : Gz));
    if (pos >= GL)
      pos -= GL;
    if (pos < 0.0)
      pos += GL;
    x[d] = pos * dx;
  }
}

static int mass_maps_write_products(FILE *file)
{
  unsigned long long pid;
  double q[3], x[3];
  int count = 0;
  for (int idx = 0; idx < MyGrids[0].total_local_size; idx++)
  {
    mass_maps_compute_from_products(idx, &pid, q, x);
    fprintf(file, "%llu %.8g %.8g %.8g %.8g %.8g %.8g\n", pid, q[0], q[1], q[2], x[0], x[1], x[2]);
    ++count;
  }
  return count;
}

/* ------------------------------------------------------------------------- */
/* write_mass_maps                                                           */
/*                                                                           */
/* Purpose:                                                                  */
/*   Entry point for MASS_MAPS particle position dump (full lattice mode).   */
/*   Primary path: write all particles using products[].                     */
/*   Fallback path: write owned fragmentation particles if products missing. */
/*                                                                           */
/* File format (per task):                                                   */
/*   ASCII, one particle per line:                                           */
/*     <ID> <q_x> <q_y> <q_z> <x_x> <x_y> <x_z>                              */
/*   Units: q_* and x_* in length units (params.InterPartDist multiples).    */
/*                                                                           */
/* Arguments:                                                                */
/*   z_start : redshift at which to output (used directly).                  */
/*   z_end   : (reserved for future interpolation / range outputs).          */
/*                                                                           */
/* Returns: 0 on success, >0 on error.                                       */
/* ------------------------------------------------------------------------- */
int write_mass_maps(double z_start, double z_end)
{
  (void)z_end;               /* currently unused */
  double redshift = z_start; /* target redshift for this dump */

  char filename[LBLENGTH];
  sprintf(filename, "massmaps.%s.z%.4f.task%d.out", params.RunFlag, redshift, ThisTask);
  FILE *file = fopen(filename, "w");
  if (!file)
  {
    printf("ERROR: Task %d could not open %s for writing\n", ThisTask, filename);
    return 1;
  }
  int nprod_written = 0;
  int fallback_frag_written = 0;
  if (mass_maps_full_available())
  {
    nprod_written = mass_maps_write_products(file);
  }
  else
  {
    fprintf(stderr, "ERROR: Task %d MASS_MAPS falling back to fragmentation view (products[] missing)\n", ThisTask);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  fclose(file);
  /* Diagnostics */
  unsigned long long written = (unsigned long long)(nprod_written ? nprod_written : fallback_frag_written);
  unsigned long long expected_local = (unsigned long long)(nprod_written ? MyGrids[0].total_local_size : subbox.Ngood);
  if (written != expected_local && internal.verbose_level >= VXX)
    printf("WARNING: Task %d wrote %llu particles (expected %llu) for mass maps (%s).\n", ThisTask, written, expected_local, (nprod_written ? "full lattice" : "fragment fallback"));

  unsigned long long global_written = 0, global_expected = 0;
  MPI_Reduce(&written, &global_written, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  unsigned long long expected_global_local = expected_local; /* each task's expected */
  MPI_Reduce(&expected_global_local, &global_expected, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (ThisTask == 0 && internal.verbose_level >= VDIAG)
    printf("[%s] MASS_MAPS: wrote positions at z=%.4f (%s, local=%llu/%llu global=%llu/%llu)\n",
           fdate(), redshift, (nprod_written ? "full lattice" : "fragment fallback"), written, expected_local, global_written, global_expected);
  return 0;
}

#endif /* MASS_MAPS */
