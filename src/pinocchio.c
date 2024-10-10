/* ######HEADER###### */

#include "pinocchio.h"
#include "def_splines.h"

void abort_code(void);
void write_cputimes(void);

int main(int argc, char **argv, char **envp)
{
  
  /* Initialize MPI */
  int got_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_level); // Hybrid MPI and OPENMP parallel
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); // gives you your rank i.e. ID
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks); // size of your pool 
  
  /* timing of the code */
  cputime.total=MPI_Wtime();
  greetings();

  /* checks that the parameter file is given in the command line */
  if (argc<2)
    {
      if (!ThisTask)
	printf("Usage: pinocchio.x parameterfile\n");
      MPI_Finalize();
      return 0;
    }

  /* initialization */
  memset(&params, 0, sizeof(param_data));
  strcpy(params.ParameterFile,argv[1]);
  if (initialization())
    abort_code();

  /******************************************/
  /*********** Standard behaviour ***********/
  /******************************************/
  if (params.ReadProductsFromDumps)
    {
      /* on request, read products from dump files and skip the first part */
      if (read_dumps())
	abort_code();
    }
  else
    {
      /* computation of collapse times and displacements */
      if (compute_fmax())
	abort_code();

      /* on request, dump products to files for skipping fmax */
      if (params.DumpProducts)
	if (dump_products())
	  abort_code();
    }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* output detailed cpu times */
  cputime.total=MPI_Wtime()-cputime.total;
  if (!ThisTask)
    write_cputimes();

  /* done */
  if (!ThisTask)
    printf("Pinocchio done!\n");
  MPI_Finalize();
  return 0;
}


void abort_code(void)
{
  printf("Task %d aborting...\n",ThisTask);
  MPI_Abort(MPI_COMM_WORLD,1);
}


void write_cputimes()
{
  printf("Total:            %14.6f\n", cputime.total);
  printf("Initialization:   %14.6f (%5.2f%%)\n", cputime.init, 100.*cputime.init/cputime.total);
  printf("  Density in PS:  %14.6f (%5.2f%%)\n", cputime.dens, 100.*cputime.dens/cputime.total);
  printf("fmax:             %14.6f (%5.2f%%)\n", cputime.fmax, 100.*cputime.fmax /cputime.total);
#ifdef TWO_LPT
  printf("  LPT:            %14.6f (%5.2f%%)\n", cputime.lpt,  100.*cputime.lpt  /cputime.total);
#endif
  printf("  Derivatives:    %14.6f (%5.2f%%)\n", cputime.deriv,  100.*cputime.deriv  /cputime.total);
  printf("    Mem transfer: %14.6f (%5.2f%%)\n", cputime.mem_transf, 100.*cputime.mem_transf  /cputime.total);
  printf("    FFTs:         %14.6f (%5.2f%%)\n", cputime.fft,  100.*cputime.fft  /cputime.total);
  printf("  Collapse times: %14.6f (%5.2f%%)\n", cputime.coll, 100.*cputime.coll /cputime.total);
  printf("    inv.collapse: %14.6f (%5.2f%%)\n", cputime.invcoll, 100.*cputime.invcoll /cputime.total);
  printf("    ellipsoid:    %14.6f (%5.2f%%)\n", cputime.ell, 100.*cputime.ell /cputime.total);
  printf("  Velocities:     %14.6f (%5.2f%%)\n", cputime.vel,  100.*cputime.vel  /cputime.total);
  printf("Total I/O:        %14.6f (%5.2f%%)\n", cputime.io,   100.*cputime.io   /cputime.total);
}
