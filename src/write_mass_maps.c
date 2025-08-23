/* MASS_MAPS minimal diagnostics skeleton (clean rewrite) */
#include "pinocchio.h"
#ifdef MASS_MAPS
#include <math.h>
#include <stdint.h>

#ifndef PLC
#error "MASS_MAPS requires PLC"
#endif

static inline double qglob_to_snapshot(double qglob, int axis)
{
  long off = subbox.start[axis] - MyGrids[0].GSstart[axis];
  double q = qglob - (double)off;
  double L = (double)MyGrids[0].GSglobal[axis];
  if (q >= L)
    q -= L;
  else if (q < 0.0)
    q += L;
  return q;
}

static void dump_core_particle_positions_z0(void)
{
  char fname[LBLENGTH];
  sprintf(fname, "pinocchio.%s.massmaps.positions.task%05d.out", params.RunFlag, ThisTask);
  FILE *fp = fopen(fname, "w");
  if (!fp)
    return;
  fprintf(fp, "# particle_id gx gy gz lagx lagy lagz dx dy dz x y z snap_x snap_y snap_z (z=0)\n");
  int istart = MyGrids[0].GSstart[_x_], jstart = MyGrids[0].GSstart[_y_], kstart = MyGrids[0].GSstart[_z_];
  int nx = MyGrids[0].GSlocal[_x_], ny = MyGrids[0].GSlocal[_y_], nz = MyGrids[0].GSlocal[_z_];
  int sx = subbox.safe[_x_], sy = subbox.safe[_y_], sz = subbox.safe[_z_];
  for (int ii = 0; ii < nx; ++ii)
  {
    int i_local = ii + sx;
    for (int jj = 0; jj < ny; ++jj)
    {
      int j_local = jj + sy;
      for (int kk = 0; kk < nz; ++kk)
      {
        int k_local = kk + sz;
        unsigned int pos = COORD_TO_INDEX(i_local, j_local, k_local, subbox.Lgwbl);
        long long gx = istart + ii, gy = jstart + jj, gz = kstart + kk;
        unsigned long long pid = 1ULL + (unsigned long long)COORD_TO_INDEX(gx, gy, gz, MyGrids[0].GSglobal);
        pos_data tp;
        set_point(i_local, j_local, k_local, (int)pos, (PRODFLOAT)1.0, &tp);
        double lag[3] = {(double)gx + SHIFT, (double)gy + SHIFT, (double)gz + SHIFT};
        double disp[3], epos[3], snap_epos[3];
        for (int a = 0; a < 3; ++a)
        {
          double qloc = q2x(a, &tp, subbox.pbc[a], (double)subbox.Lgwbl[a], ORDER_FOR_CATALOG);
          double qglob = qloc + subbox.stabl[a];
          disp[a] = qglob - lag[a];
          epos[a] = qglob * params.InterPartDist;
          double snap_q = qglob_to_snapshot(qglob, a);
          snap_epos[a] = snap_q * params.InterPartDist * (params.OutputInH100 ? params.Hubble100 : 1.0);
        }
        fprintf(fp, "%llu %lld %lld %lld %.6g %.6g %.6g %.6g %.6g %.6g %.8g %.8g %.8g %.8g %.8g %.8g\n", pid, gx, gy, gz, lag[0], lag[1], lag[2], disp[0], disp[1], disp[2], epos[0], epos[1], epos[2], snap_epos[0], snap_epos[1], snap_epos[2]);
      }
    }
  }
  fclose(fp);
}

static void summarize_core_particle_positions_z0(void)
{
  int istart = MyGrids[0].GSstart[_x_], jstart = MyGrids[0].GSstart[_y_], kstart = MyGrids[0].GSstart[_z_];
  int nx = MyGrids[0].GSlocal[_x_], ny = MyGrids[0].GSlocal[_y_], nz = MyGrids[0].GSlocal[_z_];
  int sx = subbox.safe[_x_], sy = subbox.safe[_y_], sz = subbox.safe[_z_];
  double local_sum[3] = {0}, local_sum2[3] = {0};
  double local_min[3] = {1e300, 1e300, 1e300};
  double local_max[3] = {-1e300, -1e300, -1e300};
  uint64_t local_hash = 1469598103934665603ULL;
  long long local_count = 0, local_nan = 0;
  for (int ii = 0; ii < nx; ++ii)
  {
    int i_local = ii + sx;
    for (int jj = 0; jj < ny; ++jj)
    {
      int j_local = jj + sy;
      for (int kk = 0; kk < nz; ++kk)
      {
        int k_local = kk + sz;
        unsigned int pos = COORD_TO_INDEX(i_local, j_local, k_local, subbox.Lgwbl);
        pos_data tp;
        set_point(i_local, j_local, k_local, (int)pos, (PRODFLOAT)1.0, &tp);
        double snap_epos[3];
        int bad = 0;
        for (int a = 0; a < 3; ++a)
        {
          double qglob = q2x(a, &tp, subbox.pbc[a], (double)subbox.Lgwbl[a], ORDER_FOR_CATALOG) + subbox.stabl[a];
          double snap_q = qglob_to_snapshot(qglob, a);
          double v = snap_q * params.InterPartDist * (params.OutputInH100 ? params.Hubble100 : 1.0);
          snap_epos[a] = v;
          if (!isfinite(v))
            bad = 1;
        }
        if (bad)
        {
          local_nan++;
          continue;
        }
        for (int a = 0; a < 3; ++a)
        {
          double v = snap_epos[a];
          local_sum[a] += v;
          local_sum2[a] += v * v;
          if (v < local_min[a])
            local_min[a] = v;
          if (v > local_max[a])
            local_max[a] = v;
          union
          {
            double d;
            uint64_t u;
          } bits;
          bits.d = v;
          local_hash ^= bits.u;
          local_hash *= 1099511628211ULL;
        }
        local_count++;
      }
    }
  }
  double global_sum[3], global_sum2[3], global_min[3], global_max[3];
  long long global_count = 0, global_nan = 0;
  uint64_t recv_hash = 0;
  MPI_Reduce(local_sum, global_sum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_sum2, global_sum2, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_min, global_min, 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_max, global_max, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_count, &global_count, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_nan, &global_nan, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_hash, &recv_hash, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR, 0, MPI_COMM_WORLD);
  if (ThisTask == 0)
  {
    char fname[LBLENGTH];
    sprintf(fname, "pinocchio.%s.massmaps.positions.summary.out", params.RunFlag);
    FILE *fp = fopen(fname, "w");
    if (fp)
    {
      fprintf(fp, "# Summary snapshot-frame positions z=0 run %s\n", params.RunFlag);
      fprintf(fp, "COUNT %lld\n", global_count);
      fprintf(fp, "SUM_X %.17g SUM_Y %.17g SUM_Z %.17g\n", global_sum[0], global_sum[1], global_sum[2]);
      fprintf(fp, "SUM2_X %.17g SUM2_Y %.17g SUM2_Z %.17g\n", global_sum2[0], global_sum2[1], global_sum2[2]);
      fprintf(fp, "MIN_X %.17g MIN_Y %.17g MIN_Z %.17g\n", global_min[0], global_min[1], global_min[2]);
      fprintf(fp, "MAX_X %.17g MAX_Y %.17g MAX_Z %.17g\n", global_max[0], global_max[1], global_max[2]);
      fprintf(fp, "HASH 0x%016llx\n", (unsigned long long)recv_hash);
      if (global_nan)
        fprintf(fp, "# WARNING %lld NaN\n", global_nan);
      fclose(fp);
    }
  }
}

static void summarize_core_particle_lattice(void)
{
  int istart = MyGrids[0].GSstart[_x_], jstart = MyGrids[0].GSstart[_y_], kstart = MyGrids[0].GSstart[_z_];
  int nx = MyGrids[0].GSlocal[_x_], ny = MyGrids[0].GSlocal[_y_], nz = MyGrids[0].GSlocal[_z_];
  double local_sum[3] = {0}, local_sum2[3] = {0};
  double local_min[3] = {1e300, 1e300, 1e300};
  double local_max[3] = {-1e300, -1e300, -1e300};
  uint64_t local_hash = 1469598103934665603ULL;
  long long local_count = 0;
  for (int ii = 0; ii < nx; ++ii)
    for (int jj = 0; jj < ny; ++jj)
      for (int kk = 0; kk < nz; ++kk)
      {
        double epos[3] = {(istart + ii + SHIFT) * params.InterPartDist, (jstart + jj + SHIFT) * params.InterPartDist, (kstart + kk + SHIFT) * params.InterPartDist};
        for (int a = 0; a < 3; ++a)
        {
          double v = epos[a];
          local_sum[a] += v;
          local_sum2[a] += v * v;
          if (v < local_min[a])
            local_min[a] = v;
          if (v > local_max[a])
            local_max[a] = v;
          union
          {
            double d;
            uint64_t u;
          } bits;
          bits.d = v;
          local_hash ^= bits.u;
          local_hash *= 1099511628211ULL;
        }
        local_count++;
      }
  double global_sum[3], global_sum2[3], global_min[3], global_max[3];
  long long global_count = 0;
  uint64_t recv_hash = 0;
  MPI_Reduce(local_sum, global_sum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_sum2, global_sum2, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_min, global_min, 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_max, global_max, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_count, &global_count, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_hash, &recv_hash, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR, 0, MPI_COMM_WORLD);
  if (ThisTask == 0)
  {
    char fname[LBLENGTH];
    sprintf(fname, "pinocchio.%s.massmaps.positions.lattice.summary.out", params.RunFlag);
    FILE *fp = fopen(fname, "w");
    if (fp)
    {
      fprintf(fp, "# Lattice-only core particle positions (no displacement) run %s\n", params.RunFlag);
      fprintf(fp, "COUNT %lld\n", global_count);
      fprintf(fp, "SUM_X %.17g SUM_Y %.17g SUM_Z %.17g\n", global_sum[0], global_sum[1], global_sum[2]);
      fprintf(fp, "SUM2_X %.17g SUM2_Y %.17g SUM2_Z %.17g\n", global_sum2[0], global_sum2[1], global_sum2[2]);
      fprintf(fp, "MIN_X %.17g MIN_Y %.17g MIN_Z %.17g\n", global_min[0], global_min[1], global_min[2]);
      fprintf(fp, "MAX_X %.17g MAX_Y %.17g MAX_Z %.17g\n", global_max[0], global_max[1], global_max[2]);
      fprintf(fp, "HASH 0x%016llx\n", (unsigned long long)recv_hash);
      fclose(fp);
    }
  }
}

int write_mass_maps(double z_start, double z_end)
{
  (void)z_start;
  (void)z_end;
  if (ThisTask == 0)
    printf("[%s] MASS_MAPS: diagnostics-only skeleton (positions + summaries).\n", fdate());
  dump_core_particle_positions_z0();
  summarize_core_particle_positions_z0();
  summarize_core_particle_lattice();
  return 0;
}

#endif /* MASS_MAPS */
