/* ######HEADER###### */
#include "pinocchio.h"

#define SET_WTIME
#define ASSIGN_WTIME(INIT, ACC)
#define ACCUMULATE_WTIME(INIT, ACC)

int initialize_fft(void);
int init_cosmology(void);
int generate_densities(void);
int set_plc(void);
unsigned int gcd(unsigned int, unsigned int);
int set_fft_decomposition(void);

int initialization()
{

  /* timing */
  cputime.init=MPI_Wtime();

  /* this is for gsl integration */
  workspace = gsl_integration_workspace_alloc(NWINT);
  /* this is the initialization of the random number generator */
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  /* reading of parameters from file and few other parameter initializations */
  if (set_parameters())
    return 1;

  /* call to cosmo.c: initialization of cosmological functions */
  if (initialize_cosmology())
    return 1;

  /* initialize pfft and fftw functions */
  if (initialize_fft())
    return 1;

  /* set the properties of grids and initialize FFTW quantities, including vectors */
  if (set_grids())
    return 1;

  /* computes the smoothing radii */
  if (set_smoothing())
    return 1;

  /* now it re-initializes the variance with a top-hat filter */
  WindowFunctionType=2;
  if (initialize_MassVariance())
    return 1;

  /* initializes quantities needed for the on-the-fly reconstruction of PLC */
  SET_WTIME;
  if (set_plc())
    return 1;
  ASSIGN_WTIME(partial, set_plc);

  /* computes the number of sub-boxes for fragmentation */
  SET_WTIME;
  if (set_subboxes())
    return 1;
  ASSIGN_WTIME(partial, set_subboxes);
  
  /* checks that parameters and directives are coherent */
  if (check_parameters_and_directives())
    return 1;

  /* estimates the size of output file
  if (estimate_file_size())
    return 1;
  */

  /* this barrier is set to have correct stdout in case the run cannot start */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* allocations of memory for fmax and memory tests */
  SET_WTIME;
  if (allocate_main_memory())
    return 1;
  ASSIGN_WTIME(partial, memory_allocation);

  /* initialization of fft plans */
  SET_WTIME;
  if (compute_fft_plans())
    return 1;
  ASSIGN_WTIME(partial, fft_initialization);

  /* generation of initial density field */
  if (!params.ReadProductsFromDumps) /* no generation if products are read from dumps */
    if (generate_densities())
      return 1;

  cputime.init=MPI_Wtime()-cputime.init;

  if (!ThisTask)
    {      
      dprintf(VMSG, ThisTask, "[%s] initialization done, initialization cpu time = %14.6f\n", fdate(), cputime.init);
      dprintf(VMSG, ThisTask, "\t\t set subboxes time = %14.6f\n"
	      "\t\t set plc time = %14.6f\n"
	      "\t\t memory allocation time = %14.6f\n"
	      "\t\t fft initialization time = %14.6f\n"
	      "\t\t density generation time = %14.6f\n",
	      cputime.set_subboxes, cputime.set_plc, cputime.memory_allocation,
	      cputime.fft_initialization, cputime.dens);
    }

  return 0;
}


int initialize_fft(void)
{

  /* Initialize pfft */
  pfft_init();

  /* Inititalize fftw */
  fftw_mpi_init();

  if(set_fft_decomposition())
    return 1;

  if (!ThisTask)
    dprintf(VMSG, ThisTask, "cube subdivision [%d dim]: %d x %d x %d = %d processes\n",
	    internal.tasks_subdivision_dim,
	    internal.tasks_subdivision_3D[0],
	    internal.tasks_subdivision_3D[1],
	    internal.tasks_subdivision_3D[2],
	    internal.tasks_subdivision_3D[0] *
	    internal.tasks_subdivision_3D[1] *
	    internal.tasks_subdivision_3D[2]);
  
  if ( pfft_create_procmesh(internal.tasks_subdivision_dim, MPI_COMM_WORLD, internal.tasks_subdivision_3D, &FFT_Comm) )
    {
      int all = 1;
      for(int iii = 0; iii < internal.tasks_subdivision_dim; iii++)
  	all *= internal.tasks_subdivision_3D[iii];
      
      pfft_fprintf(MPI_COMM_WORLD, stderr, "Error while creating communicator and mesh with %d processes\n", all);
      return 1;
    }

  return 0;
}


int set_parameters()
{
  int i;

  /* set default internal parameters */
  internal.verbose_level                   = VDIAG;
  internal.dump_seedplane                  = 0;
  internal.dump_kdensity                   = 0;
  internal.large_plane                     = 1;
  internal.mimic_original_seedtable        = 0;
  internal.dump_vectors                    = 0;
  internal.constrain_task_decomposition[0] = 0;
  internal.constrain_task_decomposition[1] = 0;
  internal.constrain_task_decomposition[2] = 0;
  internal.tasks_subdivision_3D[0]         = 0;
  internal.tasks_subdivision_3D[1]         = 0;
  internal.tasks_subdivision_3D[2]         = 0;

  if(read_parameter_file())
    return 1;

  /* the smallest legitimate value of MinHaloMass is 1 */
  if (params.MinHaloMass<=0)
    params.MinHaloMass=1;

  if (params.BoxInH100)
    {
      params.BoxSize_h100  = params.BoxSize;
      params.BoxSize_htrue = params.BoxSize/params.Hubble100;
    }
  else
    {
      params.BoxSize_h100  = params.BoxSize*params.Hubble100;
      params.BoxSize_htrue = params.BoxSize;
    }
  params.InterPartDist = params.BoxSize_htrue/params.GridSize[0];

  params.ParticleMass = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0 
    * pow(params.InterPartDist,3.);
  strcpy(params.DumpDir,"DumpProducts/");

  /* The Nyquist wavenumber is used in generic calls of scale-dependent growth rates */
  params.k_for_GM = PI/params.InterPartDist;

  if (!params.NumFiles)
    params.NumFiles=1;

  /* The number of files must be a divisor of the number of tasks */
  if (NTasks%params.NumFiles != 0)
    {
      while (NTasks%params.NumFiles != 0)
	params.NumFiles--;

      if (!ThisTask)
	printf("Warning: NumFiles must be a divisor of NTasks, it has been fixed to %d\n",
	       params.NumFiles);
    }


  /* inverse collapse times for the required outputs */
  for (i=0; i<outputs.n; i++)
    outputs.F[i]=1.+outputs.z[i];
  outputs.Flast=outputs.F[outputs.n-1];

  if (!ThisTask)
    {
      dprintf(VMSG, 0, "Flag for this run: %s\n\n",params.RunFlag);
      dprintf(VMSG, 0, "PARAMETER VALUES from file %s:\n",params.ParameterFile);
      dprintf(VMSG, 0, "Omega0                      %f\n",params.Omega0);
      dprintf(VMSG, 0, "OmegaLambda                 %f\n",params.OmegaLambda);    
      dprintf(VMSG, 0, "OmegaBaryon                 %f\n",params.OmegaBaryon);
      if (strcmp(params.TabulatedEoSfile,"no"))
	{
	  dprintf(VMSG, 0, "Dark Energy EoS will be read from file %s\n",params.TabulatedEoSfile);
	}
      else
	{
	  dprintf(VMSG, 0, "DE EoS parameters           %f %f\n",params.DEw0,params.DEwa);
	}

      dprintf(VMSG, 0, "Hubble100                   %f\n",params.Hubble100);
      dprintf(VMSG, 0, "Sigma8                      %f\n",params.Sigma8);
      dprintf(VMSG, 0, "PrimordialIndex             %f\n",params.PrimordialIndex);
      dprintf(VMSG, 0, "RandomSeed                  %d\n",params.RandomSeed);
      dprintf(VMSG, 0, "OutputList                  %s\n",params.OutputList);
      dprintf(VMSG, 0, "Number of outputs           %d\n",outputs.n);
      dprintf(VMSG, 0, "Output redshifts           ");
      for (i=0; i<outputs.n; i++)
	dprintf(VMSG, 0, " %f ",outputs.z[i]);
      dprintf(VMSG, 0, "\n");
      dprintf(VMSG, 0, "GridSize                    %d %d %d\n",params.GridSize[0],params.GridSize[1],params.GridSize[2]);
      dprintf(VMSG, 0, "BoxSize (true Mpc)          %f\n",params.BoxSize_htrue);
      dprintf(VMSG, 0, "BoxSize (Mpc/h)             %f\n",params.BoxSize_h100);
      dprintf(VMSG, 0, "Particle Mass (true Msun)   %g\n",params.ParticleMass);
      dprintf(VMSG, 0, "Particle Mass (Msun/h)      %g\n",params.ParticleMass*params.Hubble100);
      dprintf(VMSG, 0, "Inter-part dist (true Mpc)  %f\n",params.InterPartDist);
      dprintf(VMSG, 0, "Inter-part dist (Mpc/h)     %f\n",params.InterPartDist*params.Hubble100);
      dprintf(VMSG, 0, "MinHaloMass (particles)     %d\n",params.MinHaloMass);
      dprintf(VMSG, 0, "MinHaloMass (Msun/h)        %g\n",params.MinHaloMass*params.ParticleMass*params.Hubble100);
      dprintf(VMSG, 0, "BoundaryLayerFactor         %f\n",params.BoundaryLayerFactor);
      dprintf(VMSG, 0, "MaxMem per task (Mb)        %d\n",params.MaxMem);
      dprintf(VMSG, 0, "MaxMem per particle (b)     %f\n",params.MaxMemPerParticle);
      dprintf(VMSG, 0, "CatalogInAscii              %d\n",params.CatalogInAscii);
      dprintf(VMSG, 0, "NumFiles                    %d\n",params.NumFiles);
      dprintf(VMSG, 0, "DoNotWriteCatalogs          %d\n",params.DoNotWriteCatalogs);
      dprintf(VMSG, 0, "DoNotWriteHistories         %d\n",params.DoNotWriteHistories);
      dprintf(VMSG, 0, "WriteTimelessSnapshot       %d\n",params.WriteTimelessSnapshot);
      dprintf(VMSG, 0, "OutputInH100                %d\n",params.OutputInH100);
      dprintf(VMSG, 0, "WriteDensity                %d\n",params.WriteDensity);
      dprintf(VMSG, 0, "WriteProducts               %d\n",params.WriteProducts);
      dprintf(VMSG, 0, "DumpProducts                %d\n",params.DumpProducts);
      dprintf(VMSG, 0, "ReadProductsFromDumps       %d\n",params.ReadProductsFromDumps);
      switch(params.AnalyticMassFunction)
	{
	case 0:
	  dprintf(VMSG, 0, "Using Press & Schechter (1974) for the analytic mass function\n");
	  break;
	case 1:
	  dprintf(VMSG, 0, "Using Sheth & Tormen (2001) for the analytic mass function\n");
	  break;
	case 2:
	  dprintf(VMSG, 0, "Using Jenkins et al. (2001) for the analytic mass function\n");
	  break;
	case 3:
	  dprintf(VMSG, 0, "Using Warren et al. (2006) for the analytic mass function\n");
	  break;
	case 4:
	  dprintf(VMSG, 0, "Using Reed et al. (2007) for the analytic mass function\n");
	  break;
	case 5:
	  dprintf(VMSG, 0, "Using Crocce et al. (2010) for the analytic mass function\n");
	  break;
	case 6:
	  dprintf(VMSG, 0, "Using Tinker et al. (2008) for the analytic mass function\n");
	  break;
	case 7:
	  dprintf(VMSG, 0, "Using Courtin et al. (2010) for the analytic mass function\n");
	  break;
	case 8:
	  dprintf(VMSG, 0, "Using Angulo et al. (2012) for the analytic mass function\n");
	  break;
	case 9:
	  dprintf(VMSG, 0, "Using Watson et al. (2013) for the analytic mass function\n");
	  break;
	case 10:
	  dprintf(VMSG, 0, "Using Crocce et al. (2010) with forced universality for the analytic mass function\n");
	  break;
	default:
	  dprintf(VMSG, 0, "Unknown value for AnalyticMassFunction, Using Watson et al. (2013)\n");
	  params.AnalyticMassFunction=9;
	  break;
	}
      dprintf(VMSG, 0, "\n");

      dprintf(VMSG, 0, "\n");
      dprintf(VMSG, 0, "GENIC parameters:\n");
      dprintf(VMSG, 0, "InputSpectrum_UnitLength_in_cm %f\n",params.InputSpectrum_UnitLength_in_cm);
      dprintf(VMSG, 0, "FileWithInputSpectrum          %s\n",params.FileWithInputSpectrum);
      dprintf(VMSG, 0, "WDM_PartMass_in_kev            %f\n",params.WDM_PartMass_in_kev);
      dprintf(VMSG, 0, "\n");
    }

  /* Task 0 may have changed the value of this parameter */
  MPI_Bcast(&params.AnalyticMassFunction, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return 0;
}

int set_smoothing()
{
  int ismooth;
  double var_min, var_max, rmin;

  var_min    = pow(1.686/NSIGMA / GrowingMode(outputs.zlast,params.k_for_GM),2.0);
  rmin       = params.InterPartDist/6.;
  var_max    = MassVariance(rmin);
  Smoothing.Nsmooth = (log10(var_max)-log10(var_min))/STEP_VAR+2;

  if (Smoothing.Nsmooth<=0)
    {
      if (!ThisTask)
	dprintf(VERR, 0, "I am afraid that nothing is predicted to collapse in this configuration.\nI will work with no smoothing\n");
      Smoothing.Nsmooth=1;
    }

  if (!ThisTask)
    {
      printf("\nSMOOTHING RADII\n");
      printf("Min variance: %f12.6, max variance: %f12.6, number of smoothing radii: %d\n",
	     var_min,var_max,Smoothing.Nsmooth);
    }
  Smoothing.Radius      =(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.Variance    =(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.TrueVariance=(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  if (Smoothing.Radius==0x0 || Smoothing.Variance==0x0 || Smoothing.TrueVariance==0x0)
    {
      printf("ERROR on task %d: allocation of Smoothing failed\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  for (ismooth=0; ismooth<Smoothing.Nsmooth-1; ismooth++)
    {
      Smoothing.Variance[ismooth] = pow(10., log10(var_min)+STEP_VAR*ismooth);
      Smoothing.Radius[ismooth]   = Radius(Smoothing.Variance[ismooth]);
    }
  Smoothing.Radius[ismooth]   = 0.0;
  Smoothing.Variance[ismooth] = var_max;

  if (!ThisTask)
    for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
      printf("           %2d)  Radius=%10f, Variance=%10f\n",ismooth+1,Smoothing.Radius[ismooth],Smoothing.Variance[ismooth]);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}


int generate_densities()
{

  cputime.dens=MPI_Wtime();

  if (!ThisTask)
    dprintf(VMSG, 0, "[%s] Generating density in Fourier space\n",fdate());

  int igrid;
  for (igrid=0; igrid<Ngrids; igrid++)
    if (GenIC_large(igrid))
      return 1;

  cputime.dens = MPI_Wtime()-cputime.dens;
    if (!ThisTask)
      dprintf(VMSG, 0, "[%s] Done generating density in Fourier space, cputime = %f s\n",fdate(), cputime.dens);

  return 0;
}


int set_grids()
{
  /* initialization of fftw quantities on grids (one for the moment) */

  int igrid, dim;

  Ngrids=1;

  MyGrids=(grid_data*)malloc(Ngrids * sizeof(grid_data));

  for (dim=0; dim<3; dim++)
    MyGrids[0].GSglobal[dim] = params.GridSize[dim];
  
  MyGrids[0].Ntotal = (unsigned long long)MyGrids[0].GSglobal[_x_] * 
    (unsigned long long)MyGrids[0].GSglobal[_y_] * 
    (unsigned long long)MyGrids[0].GSglobal[_z_];

  MyGrids[0].BoxSize = params.BoxSize_htrue;
  MyGrids[0].lower_k_cutoff=0.;
  MyGrids[0].upper_k_cutoff=NYQUIST * PI;

  /* allocates pointers */
  cvector_fft=(pfft_complex**)malloc(Ngrids * sizeof(fftw_complex*));
  rvector_fft=(double**)malloc(Ngrids * sizeof(double*));

  kdensity=(double**)malloc(Ngrids * sizeof(double*));
  density=(double**)malloc(Ngrids * sizeof(double*));
  first_derivatives=(double***)malloc(Ngrids * sizeof(double**));
  second_derivatives=(double***)malloc(Ngrids * sizeof(double**));
  VEL_for_displ=(double**)malloc(3 * sizeof(double*));
#ifdef TWO_LPT
  VEL2_for_displ=(double**)malloc(3 * sizeof(double*));
#endif
  for (igrid=0; igrid<Ngrids; igrid++)
    {
      first_derivatives[igrid]=(double**)malloc(3 * sizeof(double*));
      second_derivatives[igrid]=(double**)malloc(6 * sizeof(double*));
    }
  /* moved to GenIC */
  /* seedtable=(unsigned int**)malloc(Ngrids * sizeof(unsigned int*)); */
 
  for (igrid=0; igrid<Ngrids; igrid++)
    if (set_one_grid(igrid))
      return 1;

  /* Task 0 broadcasts its number of fft particles, that becomes the reference */
  unsigned int PPT;
  if (!ThisTask)
    PPT=MyGrids[0].total_local_size;

  MPI_Bcast(&PPT, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MyGrids[0].ParticlesPerTask=PPT;

  if ((int)(PPT * params.MaxMemPerParticle / MBYTE + 1.0) > params.MaxMem)
    {
      if (!ThisTask)
	{
	  printf("ERROR: MaxMem of %d Mb per task is insufficient to store %d bytes per particle\n",
		 params.MaxMem, (int)params.MaxMemPerParticle);
	  printf("       please increase MaxMem to at least %d\n",(int)(PPT * params.MaxMemPerParticle / MBYTE + 1.0));
	}
      return 1;
    }

  return 0;
}

int set_plc()
{
  if (!ThisTask)
    printf("PLC flag at compilation was not set, no Past Light Cone output will be given\n\n");

  return 0;
}

/* division in sub-boxes */
int set_subboxes()
{

  int i,j,k, i1=0,j1=0,k1=0, NN,NN1,N1,N2,N3;
  unsigned long long int surface,this,tt;
  double size,sizeG,cc;

  /* mass of the largest halo expected in the box */
  params.Largest=1.e18;
  cc=1./pow(params.BoxSize_htrue,3.0);
  double aa=AnalyticMassFunction(params.Largest,outputs.zlast);
  while (aa*params.Largest<cc)
    {
      params.Largest*=0.99;
      aa=AnalyticMassFunction(params.Largest,outputs.zlast);
    }
  size=SizeForMass(params.Largest);
  sizeG=size/params.InterPartDist;
  
  /*  
      The number of loadable subbox particles is equal 
      to the number allowed by the specified MaxMemPerParticle.
      The boundary layer is set to its maximum value.
  */

  /* finds the optimal number of sub-boxes to use for the fragmentation */


  surface=MyGrids[0].Ntotal;
  for (k=1; k<=NTasks; k++)
    for (j=1; j<=NTasks/k; j++)
      for (i=1; i<=NTasks/k/j; i++)
	/* the three indices must be exact divisors of the three grid lengths */
	if (i*j*k==NTasks)
	  {
	    /* number of particles in the sub-box */
	    N1 = find_length(MyGrids[0].GSglobal[_x_],i,0);
	    N2 = find_length(MyGrids[0].GSglobal[_y_],j,0);
	    N3 = find_length(MyGrids[0].GSglobal[_z_],k,0);

	    this = (unsigned long long int)(i>1? 2*(N2*N3) : 0) + 
	      (unsigned long long int)(j>1? 2*(N1*N3) : 0) +
	      (unsigned long long int)(k>1? 2*(N1*N2) : 0);
	    tt=this;
	    if (N1/2 < sizeG)
	      this+=(unsigned long long int)((double)tt*pow(2*sizeG/(double)N1,2.0));
	    if (N2/2 < sizeG)
	      this+=(unsigned long long int)((double)tt*pow(2*sizeG/(double)N2,2.0));
	    if (N3/2 < sizeG)
	      this+=(unsigned long long int)((double)tt*pow(2*sizeG/(double)N3,2.0));

	    if (this<surface)
	      {
		surface=this;
		i1=i; 
		j1=j; 
		k1=k; 
	      }
	  }

  subbox.nbox[_x_]=i1;
  subbox.nbox[_y_]=j1;
  subbox.nbox[_z_]=k1;

  /* mybox is the box assigned to the task */
  NN=subbox.nbox[_y_]*subbox.nbox[_z_];
  if (NN==0)
    {
      printf("ERROR: I could not find a valid subbox subdivision\n");
      printf("       subbox.nbox = [%d,%d,%d]\n",i1,j1,k1);
      printf("       please try again with a different number of tasks\n");
      return 1;
    }

  subbox.mybox[_x_]=ThisTask/NN;
  NN1=ThisTask-subbox.mybox[_x_]*NN;
  subbox.mybox[_y_]=NN1/subbox.nbox[_z_];
  subbox.mybox[_z_]=NN1-subbox.mybox[_y_]*subbox.nbox[_z_];

  subbox.Lgrid[_x_] = find_length(MyGrids[0].GSglobal[_x_],subbox.nbox[_x_],subbox.mybox[_x_]);
  subbox.Lgrid[_y_] = find_length(MyGrids[0].GSglobal[_y_],subbox.nbox[_y_],subbox.mybox[_y_]);
  subbox.Lgrid[_z_] = find_length(MyGrids[0].GSglobal[_z_],subbox.nbox[_z_],subbox.mybox[_z_]);

  subbox.pbc[_x_] = (subbox.nbox[_x_]==1);
  subbox.pbc[_y_] = (subbox.nbox[_y_]==1);
  subbox.pbc[_z_] = (subbox.nbox[_z_]==1);

  int BB = (int)(params.BoundaryLayerFactor*sizeG+1);
  subbox.safe[_x_] = (subbox.pbc[_x_] ? 0 : (BB > MyGrids[0].GSglobal[_x_]/2 - subbox.Lgrid[_x_]/2 - 1 ? MyGrids[0].GSglobal[_x_]/2 - subbox.Lgrid[_x_]/2 - 1 : BB));
  subbox.safe[_y_] = (subbox.pbc[_y_] ? 0 : (BB > MyGrids[0].GSglobal[_y_]/2 - subbox.Lgrid[_y_]/2 - 1 ? MyGrids[0].GSglobal[_y_]/2 - subbox.Lgrid[_y_]/2 - 1 : BB));
  subbox.safe[_z_] = (subbox.pbc[_z_] ? 0 : (BB > MyGrids[0].GSglobal[_z_]/2 - subbox.Lgrid[_z_]/2 - 1 ? MyGrids[0].GSglobal[_z_]/2 - subbox.Lgrid[_z_]/2 - 1 : BB));

  subbox.Lgwbl[_x_] = subbox.Lgrid[_x_] + 2*subbox.safe[_x_];
  subbox.Lgwbl[_y_] = subbox.Lgrid[_y_] + 2*subbox.safe[_y_];
  subbox.Lgwbl[_z_] = subbox.Lgrid[_z_] + 2*subbox.safe[_z_];
  unsigned long long MySize = (long long)subbox.Lgwbl[_x_] * (long long)subbox.Lgwbl[_y_] * (long long)subbox.Lgwbl[_z_];
  while (MySize > (unsigned long long)1<<31)
    {
      subbox.safe[_x_] -=1;
      subbox.safe[_y_] -=1;
      subbox.safe[_z_] -=1;
      subbox.Lgwbl[_x_] = subbox.Lgrid[_x_] + 2*subbox.safe[_x_];
      subbox.Lgwbl[_y_] = subbox.Lgrid[_y_] + 2*subbox.safe[_y_];
      subbox.Lgwbl[_z_] = subbox.Lgrid[_z_] + 2*subbox.safe[_z_];
      MySize = (unsigned long long)subbox.Lgwbl[_x_] * (unsigned long long)subbox.Lgwbl[_y_] * (unsigned long long)subbox.Lgwbl[_z_];
    }
  
  subbox.start[_x_] = find_start(MyGrids[0].GSglobal[_x_],subbox.nbox[_x_],subbox.mybox[_x_]);
  subbox.start[_y_] = find_start(MyGrids[0].GSglobal[_y_],subbox.nbox[_y_],subbox.mybox[_y_]);
  subbox.start[_z_] = find_start(MyGrids[0].GSglobal[_z_],subbox.nbox[_z_],subbox.mybox[_z_]);

  subbox.stabl[_x_] = subbox.start[_x_] - subbox.safe[_x_];
  subbox.stabl[_y_] = subbox.start[_y_] - subbox.safe[_y_];
  subbox.stabl[_z_] = subbox.start[_z_] - subbox.safe[_z_];

  /* 
     Npart: total number of particles in the whole sub-volume 
     Ngood: total number of particles in the well reconstructed region
     Npredpeaks: a guess of the maximum number of peaks in the subbox
     Nalloc: number of particles for which memory has been allocated (set in organize_main_memory)
     Nstored: number of actually stored particles
  */

  subbox.Npart = subbox.Lgwbl[_x_] * subbox.Lgwbl[_y_] * subbox.Lgwbl[_z_];
  subbox.Ngood = subbox.Lgrid[_x_] * subbox.Lgrid[_y_] * subbox.Lgrid[_z_];
  /* this is a prediction of the number of peaks that will be found */
  subbox.PredNpeaks = (int)(MyGrids[0].ParticlesPerTask/6 * params.PredPeakFactor);
  subbox.Nstored = 0;
  /* this is the size of frag_map*/
  subbox.maplength = subbox.Npart/UINTLEN + (subbox.Npart%UINTLEN!=0);
  if ( (subbox.Nalloc = organize_main_memory()) == 0 )
    {
      fflush(stdout);
      if (!ThisTask)
	printf("organize_main_memory returned an invalid Nalloc, exiting\n");
      return 1;
    }

  /* messagges */
  if (!ThisTask)
    {
      printf("\n");
      printf("FRAGMENTATION:\n");
      printf("Reference number of particles:         %d\n",MyGrids[0].ParticlesPerTask);
      printf("Requested bytes per particle:          %d\n",(int)params.MaxMemPerParticle);
      printf("Number of sub-boxes per dimension:     %d %d %d\n",subbox.nbox[_x_],subbox.nbox[_y_],subbox.nbox[_z_]);
      printf("Periodic boundary conditions:          %d %d %d\n",subbox.pbc[_x_],subbox.pbc[_y_],subbox.pbc[_z_]);
      printf("Core 0 will work on a grid:            %d %d %d\n",subbox.Lgwbl[_x_],subbox.Lgwbl[_y_],subbox.Lgwbl[_z_]);
      printf("The resolved box will be:              %d %d %d\n",subbox.Lgrid[_x_],subbox.Lgrid[_y_],subbox.Lgrid[_z_]);
      printf("Boundary layer:                        %d %d %d\n",subbox.safe[_x_],subbox.safe[_y_],subbox.safe[_z_]);
      printf("Boundary layer factor:                 %f\n",params.BoundaryLayerFactor);
      printf("Number of total particles for core 0:  %d\n",subbox.Npart);
      printf("Number of good particles for core 0:   %d\n",subbox.Ngood);
      printf("Particles that core 0 will allocate:   %d\n",subbox.Nalloc);
      printf("Allowed overhead for boundary layer:   %f\n",(float)subbox.Nalloc/(float)MyGrids[0].ParticlesPerTask);
      printf("Largest halo expected in this box at z=%f: %e Msun\n",
  	     outputs.zlast, params.Largest);
      printf("   its Lagrangian size: %f Mpc (%6.2f grid points)\n",size,sizeG);
      printf("   this requires a boundary layer of %6.2f grid points \n",sizeG*params.BoundaryLayerFactor);

#ifndef BL_GRANDISSIMA 
      if ((!subbox.pbc[_x_] && params.BoundaryLayerFactor*sizeG>subbox.safe[_x_]) || 
	  (!subbox.pbc[_y_] && params.BoundaryLayerFactor*sizeG>subbox.safe[_y_]) || 
	  (!subbox.pbc[_z_] && params.BoundaryLayerFactor*sizeG>subbox.safe[_z_]))
	{
	  printf("WARNING: the boundary layer on some dimension is smaller than the predicted size of the largest halos\n");
	  printf("         times the BoundaryLayerFactor, the most massive halos may be inaccurate\n");
	}
#endif
    }

  /* initialization of quantities required by compute_mf */
   if (params.OutputInH100)
    mf.hfactor=params.Hubble100;
  else
    mf.hfactor=1.0;
  mf.hfactor4=pow(mf.hfactor,4.);
  mf.vol=(double)MyGrids[0].Ntotal*pow(params.InterPartDist,3.0);
  mf.mmin=log10(params.MinHaloMass*params.ParticleMass)-0.001*DELTAM;
  mf.mmax=log10(params.Largest)+3.0*DELTAM;
  mf.NBIN = (int)((mf.mmax-mf.mmin)/DELTAM) +1;
  mf.ninbin=(int*)calloc(mf.NBIN,sizeof(int));
  mf.ninbin_local=(int*)calloc(mf.NBIN,sizeof(int));
  mf.massinbin=(double*)calloc(mf.NBIN,sizeof(double));
  mf.massinbin_local=(double*)calloc(mf.NBIN,sizeof(double));

  /* messages */
  if (!ThisTask)
    {
      printf("\nThe mass function will be computed from Log M=%f to Log M=%f (%d bins)\n",
      	     mf.mmin, mf.mmax, mf.NBIN);
      printf("\n");
      fflush(stdout);      
    }

  return 0;
}


int find_start(int L,int n,int ibox)
{
  int LL,MM;

  if (n==1)
    return 0;
  else
    {
      LL=L/n;
      MM=L%n;
      if (ibox==0)
        return 0;
      else if (ibox<=MM)
        return ibox*(LL+1);
      else
        return ibox*LL+MM;
    }

}


int find_length(int L, int n, int ibox)
{
  /* finds the length of a subbox, given the grid length L, 
     the number of subboxes n and the subbox id ibox */

  int LL,MM;

  if (n==1) 
    return L;
  else
    {
      LL=L/n;
      MM=L%n;
      if (ibox<MM)
        return LL+1;
      else
	return LL;
    }
}
unsigned int gcd(unsigned int u, unsigned int v)
// this version of greatest common divisor taken
// from Daniel Lemire's blog
// lemire.me/blog/2013/12/26/fastest-way-to-compute-the-greatest-common-divisor/
{
    if (u == 0) return v;
    if (v == 0) return u;
    int shift = __builtin_ctz(u | v);
    u   >>= __builtin_ctz( u );
    do {
        v >>= __builtin_ctz( v );
        if (u > v) {
            unsigned int t = v;
            v = u;
            u = t;
        }  
        v = v - u;
    } while (v != 0);
    return u << shift;
}


int set_fft_decomposition(void)
{

  /* initialize task mesh for pfft */
  /* it's up to you to decide HOW to subdivide work in 3D, and then to store it in task_subdivision_3D*/

  int decomposition_done = 0;
  
  if(  internal.constrain_task_decomposition[0] +
       internal.constrain_task_decomposition[1] +
       internal.constrain_task_decomposition[2] > 0)

    // some constraints about how to decompose fft are set in parameter file
    {
      // --- check trivial errors
      // just to be sure about trivial typos, check that none is < 0
      if (  internal.constrain_task_decomposition[0] < 0 ||
	    internal.constrain_task_decomposition[1] < 0 ||
	    internal.constrain_task_decomposition[2] < 0 )
	{
	  dprintf(VXERR, 0, "you can't constraint FFt decomposition with negative values\n");
	  return 1;
	}
      
      if ( internal.constrain_task_decomposition[0] == 0 )
	{	
	  dprintf(VXERR, 0, "you can't constraint FFt decomposition leaving first dimension to 0\n");
	  return 1;
	}
      // -------------------------
      
      // set first dimension
      internal.tasks_subdivision_3D[0] = internal.constrain_task_decomposition[0];
      internal.tasks_subdivision_3D[1] = internal.tasks_subdivision_3D[2] = 1;
      decomposition_done = 1;
      
      if( internal.constrain_task_decomposition[0] == NTasks )
	{
	  // all tasks are in dimension 1
	  internal.tasks_subdivision_dim = 1;
	  decomposition_done = 3;
	}
      
      else if( internal.constrain_task_decomposition[1] > 0 )
	{
	  // --- check trivial errors
	  if(internal.constrain_task_decomposition[0]*internal.constrain_task_decomposition[1] > NTasks )
	    {
	      dprintf(VXERR, 0, "you specified a wrong fft decomposition: Dim0 x Dim1 = %d > %d tasks\n",
		      internal.constrain_task_decomposition[0]*internal.constrain_task_decomposition[1], NTasks);
	      return 1;
	    }
	  // ------------------------
	  
	  internal.tasks_subdivision_3D[1] = internal.constrain_task_decomposition[1];
	  
	  if(internal.constrain_task_decomposition[0]*internal.constrain_task_decomposition[1] == NTasks)
	    {
	      // all tasks are in dimension 1 and 2
	      internal.tasks_subdivision_dim = 2;
	      internal.tasks_subdivision_3D[2] = 1;
	      decomposition_done = 3;
	    }
	  else
	    {
	      internal.tasks_subdivision_dim = 3;
	      decomposition_done = 3;
	      
	      internal.tasks_subdivision_3D[2] = NTasks /(internal.tasks_subdivision_3D[0] * internal.tasks_subdivision_3D[1]);
	      
	      if(internal.constrain_task_decomposition[2] > 0 &&
		 internal.tasks_subdivision_3D[2] != internal.constrain_task_decomposition[2])
		{
		  dprintf(VXERR, 0, "you specified a wrong fft decomposition: dim2 should be %d instead of %d\n", internal.tasks_subdivision_3D[2], internal.constrain_task_decomposition[2]);
		  return 1;
		}
	    }
	}
      else // constrain_task_decomposition[1] > 0
	decomposition_done = 1;
      
      // close if constrain_task_decomposition[1] > 0
      
    } // close constrain_task_decomposition initial if
  

  if(decomposition_done < 3)
    {
      // decomposition is still to be made or completed
      // we try to use as less dimensions as possible, maximizing contiguity

      
      if(decomposition_done == 1)
	// only the first dimension has been specified in the param file,
	// but with a number of tasks smaller than NTasks
	{	  
	  internal.tasks_subdivision_3D[1] = NTasks / internal.tasks_subdivision_3D[0];
	  internal.tasks_subdivision_3D[2] = 1;
	  internal.tasks_subdivision_dim = 2;
	  
	  if(NTasks % internal.tasks_subdivision_3D[0])
	    {
	      dprintf(VXERR, 0, "you specified a wrong fft decomposition\n");
	      return 1;
	    }  	  
	}
      else
	// no dimension has been constrained in the param file
	{
	  // NOTE : no non-cubic grids, no multi-grids
	  
	  int Ngrid = params.GridSize[0];

	  if( NTasks <= Ngrid )
	    // prefer 1D decomposition for the most obvious case
	    // that minimize communications in FFTs
	    {
	      internal.tasks_subdivision_dim = 1;
	      internal.tasks_subdivision_3D[0] = NTasks;
	      internal.tasks_subdivision_3D[2] = internal.tasks_subdivision_3D[1] = 1;
	      return 0;
	    }

	  int Ngrid2 = Ngrid * Ngrid / (DECOMPOSITION_LIMIT_FACTOR_2D * DECOMPOSITION_LIMIT_FACTOR_2D);

	  if( NTasks <= Ngrid2)
	    // check whether exact 2D pencil decomposition
	    {
	      unsigned GCD   = gcd( Ngrid, NTasks );
	      unsigned GCD_2 = gcd( Ngrid, (NTasks / GCD) );

	      if( GCD * GCD_2 != NTasks)
		// no exact decomposition is possible,
		// revert to 1D decomposition
		{
		  internal.tasks_subdivision_dim = 1;
		  internal.tasks_subdivision_3D[0] = NTasks;
		  internal.tasks_subdivision_3D[2] = internal.tasks_subdivision_3D[1] = 1;
		  return 0;
		}
	      else
		{
		  internal.tasks_subdivision_dim = 2;
		  internal.tasks_subdivision_3D[0] = GCD;
		  internal.tasks_subdivision_3D[1] = GCD_2;
		  return 0;
		}

	    }  // close if( NTasks < Ngrid2)

	  else
	    // try 3d decomposition
	    {
	      int Ngrid_limit = Ngrid / DECOMPOSITION_LIMIT_FACTOR_2D;
	      
	      unsigned GCD   = gcd( Ngrid_limit, NTasks );
	      unsigned GCD_2 = gcd( Ngrid_limit, (NTasks / GCD) );
	      
	      if( NTasks % (GCD * GCD_2) )
		{
		  dprintf(VXERR, 0, "3D decomposition is not possible\n");
		  return 1;
		}

	      internal.tasks_subdivision_dim = 3;
	      internal.tasks_subdivision_3D[0] = GCD;
	      internal.tasks_subdivision_3D[1] = GCD_2;

	      internal.tasks_subdivision_3D[2] = NTasks / (GCD * GCD_2);
	      return 0;
	    }
	  
	}
    }

  return 0;
}


int check_parameters_and_directives(void)
{
  if (params.WriteTimelessSnapshot || params.WriteDensity)
    {
      if (!ThisTask)
	printf("ERROR: to produce a snapshot you have to compile with SNAPSHOT (DEPRECATED) directive\n");
      return 1;
    }
  return 0;
}


void greetings(void)
{
  /* This is a list of messages to declare the most relevant precompiler directives in the stdout */

  if (!ThisTask)
    {
      printf("[%s] This is pinocchio V5.0, running on %d MPI tasks\n\n",fdate(),NTasks);

#ifdef TWO_LPT
#ifndef THREE_LPT
      printf("This version uses 2LPT displacements\n");
#else
      printf("This version uses 3LPT displacements\n");
#endif
#else
      printf("This version uses Zeldovich displacements\n");
#endif
#ifdef NORADIATION
      printf("Radiation is not included in the Friedmann equations\n");
#else
      printf("Radiation is included in the Friedmann equations\n");
#endif
#ifdef LIGHT_OUTPUT
      printf("Catalogs will be written in the light version\n");
#endif

      printf("Ellipsoidal collapse will be computed as Monaco (1995)\n");

#ifdef NO_RANDOM_MODULES
      printf("Initial conditions will be generated with non-random modules of the Fourier modes\n");
#endif

      printf("\n");

    }
}
