/*****************************************************************
 *                        PINOCCHI0  V5.0                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco
 Copyright (C) 2016
 
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include "pinocchio.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <immintrin.h>
#include <assert.h>

/*------------------------------------------------------- Macros declaration --------------------------------------------------------*/


#define SMALL ((double)1.e-20)
#define INV_3 (1.0 / 3.0)

/*------------------------------------------------------- Functions definition -------------------------------------------------------*/



/* Ellipsoidal solution at 3rd perturbative order in two different ways */

double ell_classic(const  int, const double, const double, const double);  
__attribute__((always_inline)) double ell_sng                   ( int, double, double, double);  

/* Collapase time calculation */
double ell(const int, const double, const double, const double); 

/* Calculation of inverse collpase time */

double inverse_collapse_time(const int,
			     const double *const restrict,
			     double *const restrict,
			     double *const restrict,
			     double *const restrict,
			     int *const restrict);

/* Order inverse collpase time */

void ord(double *const restrict, double *const restrict, double *const restrict);

/*------------------------------------------------------- Global Variables declaration ----------------------------------------------------*/


/* Declaring the replicated variables that will be use by each single thread */
/* NOTE: Threadprivate directive specifies that variables are replicated, with each thread having its own copy. It's a declarative directive */

#define cputime_ell            cputime.ell

/* Declaring smoothing radius */

int ismooth;

/*------------------------------------------------------- Functions implementation --------------------------------------------------------*/

/* Classical ellipsoidal collapse solution */

double ell_classic(const int    ismooth,
		   const double l1,
		   const double l2,
		   const double l3)
{
  /*
    This routine computes the smallest non-negative solution of the 3rd
    order equation for the ellipsoid, and corrects it to reproduce the
    spherical collapse correctly.
  */

  /* Local variables declaration */
  double ell;
  const double del = (l1 + l2 + l3);
  const double det = (l1 * l2 * l3);    

  /* Vanishing lambda1 eigenvalue case */
  if (fabs(l1) < SMALL)
    {
      ell = -0.1;
    }
  else  /* Not vanishing lambda1 eigenvalue case */
    {
      const double den = det / 126. + 5. * l1 * del * (del - l1) / 84.;
      /* Check 1st perturbative order conditions */
      if (fabs(den) < SMALL)
	{
	  if (fabs(del - l1) < SMALL)
	    {
	      ell = ((l1 > 0.0) ? (1.0 / l1) : -0.1); /* Zel'dovich approximation */
            }
	  else
	    {
	      /* Check 2nd perturbative order conditions */
	      const double dis = (7.0 * l1 * (l1 + 6.0 * del));
	      ell = ((dis < 0.0) ? -0.1 : (7. * l1 - sqrt(dis)) / (3. * l1 * (l1 - del)));
	      ell = ((ell < 0.0) ? -0.1 : ell);
            }
        } /* 1st perturbative order conditions */
      else
	{
	  /* 3rd order perturbative solution. For more details about the equations implemented, see Monaco 1996a */
            
	  /* Intermediate values */
	  const double rden = (1.0 / den);
	  const double a1 = 3. * l1 * (del - l1) / 14. * rden;
	  const double a1_2 = a1 * a1;
	  const double a2 = l1 * rden;
	  const double a3 = -1.0 * rden;

	  /* The collapse time b_c will be a combination of R, Q, and D == den */
	  const double q = (a1_2 - 3. * a2) / 9.;
	  const double r = (2. * a1_2 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
	  const double r_2_q_3 = r * r - q * q * q;

	  /* Check 3rd perturbative order conditions */

	  /* ---------------- Case 1 --------------- */
	  /* If R^2 - Q^2 > 0, which is valid for spherical and quasi-spherical perturbations */
            
	  /* 3rd order solution */
	  if (r_2_q_3 > 0)
	    {
	      const double fabs_r = fabs(r);
	      const double sq = pow(sqrt(r_2_q_3) + fabs_r, 0.333333333333333);
	      ell = -fabs_r / r * (sq + q / sq) - a1 / 3.;
	      ell = ((ell < 0.0) ? -0.1 : ell);
            }

	  /* ---------------- Case 2 --------------- */
	  /* The solution has to be chosen as the smallest non-negative one between s1, s2, and s3 */
            
	  /* 3rd order solution */
	  else
	    {
	      const double sq = (2.0 * sqrt(q));
	      const double t = acos(2.0 * r / q / sq);
	      const double a1_inv_3 = (a1 * INV_3);

	      double s1 = -sq * cos(t * INV_3) - a1_inv_3;
	      s1 = ((s1 < 0.0) ? 1.e10 : s1);

	      double s2 = -sq * cos((t + 2. * PI) * INV_3) - a1_inv_3;
	      s2 = ((s2 < 0.0) ? 1.e10 : s2);

	      double s3 = -sq * cos((t + 4. * PI) * INV_3) - a1_inv_3;
	      s3 = ((s3 < 0.0) ? 1.e10 : s3);

	      ell = (s1 < s2  ? s1 : s2);
	      ell = (s3 < ell ? s3 : ell);
	      ell = ((ell == 1.e10) ? -0.1 : ell);
            }
        } /* 3rd order perturbative solution */
    } /* Not vanishing lambda1 eigenvalue case */
    
  if ((del > 0.) && (ell > 0.))
    {
      const double inv_del = 1.0 / del;
      ell += -.364 * inv_del * exp(-6.5 * (l1 - l2) * inv_del - 2.8 * (l2 - l3) * inv_del);
    }

  return ell;
}

/* Ellipsoidal collapse following Nadkarni-Ghosh & Singhal (2016) */

/* We follow the dynamics of triaxial collapse in terms of eigenvalues of the deformation tensor (lambda_a), the velocity derivative tensor (lambda_v) and the gravity Hessian 
(lambda_d). The idea is: starting from BM96, where the dynamic is characterized by the evolution of the three principal axes (the physical coordinates is r=a_i(t)q(i = 1,2,3)). 
In this case the evolution of the ellipse is completely determined once six parameters are known: the three axes lengths and their velocities at some initial epoch a_init. 
An alternate description of the ellipse can be given by a set of nine dimensionless parameters lambda_a, lambda_v, lambda_d. In this description when an axis is collpasing 
lamda_a ----> 1, whereas for an expanding axes lambda_a -----> -inf.
delta = lambda_d_1 + lambda_d_2 + lambda_d_3 
The nine (dimensionless) eigenvalues completely characterize the density, velocity and shape perturbations in this model of ellipsoidal collapse. 
For the net system for the nine eigenvalues see Nadkarni-Ghosh & Singhal (2016) pag. 5 */


/* We have then to solve a system of differential equations          */
/* System of ODEs: specify f_i(t) = r.h.s. of differential equations */
/* f[i] = d(l_a)/da; f[i+3] = d(l_v)/da; d(l_d)/da                   */

int sng_system(double t, const double y[], double f[], void *sng_par) {

    /* Local variables declaration */
	int i,j;
	double sum;

	/* Needed cosmological parameters calculation at the given z */
	double omegam = OmegaMatter(1. / t-1.);  /* Cosmological mass density parameter as a function of redshift DIMENSIONLESS. From cosmo.c */
	double omegal = OmegaLambda(1. / t-1.);  /* Cosmological mass density parameter as a function of redshift DIMENSIONLESS. From cosmo.c */

	/* y array contains the eigenvalues lambda_a, lambda_v and lambda_d */
	/* In this case y[6] = lambda_d_1, y[7] = lambda_d_2, y[8] = lambda_d_3 */
	double delta = y[6] + y[7] + y[8];

    /* Calculating r.h.s sums of lambda_d_i evolution */
	/* Loop from 0 -----> 3 beacause lambda_i has 3 component (it's the same for lambda_a and lambda_v) */
	
	/* ---------------- Sum  i != j --------------- */
	 for (i = 0; i < 3; i++) {
        sum = 0.;
        for (j = 0; j < 3; j++) {
            if (i == j || y[i] == y[j]) {
                continue;
            } else {
                sum += (y[j + 6] - y[i + 6]) * ((1. - y[i]) * (1. - y[i]) * (1. + y[i + 3]) -
                       (1. - y[j]) * (1. - y[j]) * (1. + y[j + 3])) /
                       ((1. - y[i]) * (1. - y[i]) - (1. - y[j]) * (1. - y[j]));
            }
        }

		/* ---------------- Gathering r.h.s of lambda_a_i equation --------------- */
        f[i] = (y[i + 3] * (y[i] - 1.0)) / t;

		  /* ---------------- Gathering r.h.s of lambda_v_i equation --------------- */
        /* NOTE: here is the part where models of modified gravity can be introduced */
        f[i + 3] = (0.5 * (y[i + 3] * (omegam - 2.0 * omegal - 2.0) 
	
                       - 3.0 * omegam * y[i + 6]
                       - 2.0 * y[i + 3] * y[i + 3])) / t;

		/* ---------------- Gathering r.h.s of lambda_d_i equation --------------- */
        f[i + 6] = ((5. / 6. + y[i + 6]) *
                    ((3. + y[3] + y[4] + y[5]) - (1. + delta) / (2.5 + delta) * (y[3] + y[4] + y[5])) -
                    (2.5 + delta) * (1. + y[i + 3]) + sum) / t;
    }

    return GSL_SUCCESS;
}


/* Solving the ODE system for the ellipsoidal collapse following SNG */

inline double  ell_sng(int ismooth, double l1, double l2, double l3) {

	/* Needed variables for the integrations step */
	double ode_param, hh = 1.e-6;  // hh = initial step-size, ode_param = arbitrary parameters of the system == Smoothing_radius in our case

    /* Setting integration time */
	double amin = 1.e-5, amax = 5.0;
	double mya = amin;
    
	/* Step type: Explicit embedded Runge-Kutta-Fehlberg (4, 5) method */
	const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

	/* Newly allocated instance of a stepping function of type T for a system of dim = 9 */
	gsl_odeiv2_step    *ode_s     = gsl_odeiv2_step_alloc(T,9); 

	/* The control function examines the proposed change to the solution produced by a stepping function and attempts to determine the optimal step-size for a user-specified level of error.
	The standard control object is a four parameter heuristic based on absolute and relative errors eps_abs and eps_rel,and scaling factors a_y and a_dydt 
	for the system state y(t) and derivatives y'(t) respectively. */
	gsl_odeiv2_control *ode_c     = gsl_odeiv2_control_standard_new(1.0e-6, 1.0e-6, 1.0, 1.0);

	/* The evolution function combines the results of a stepping function and control function to reliably advance the solution forward one step using an acceptable step-size.
	This function returns a pointer to a newly allocated instance of an evolution function for a system of dim = 9 */
	gsl_odeiv2_evolve  *ode_e     = gsl_odeiv2_evolve_alloc(9);

	/* A system of equations is defined using the gsl_odeiv2_system datatype 
	This data type defines a general ODE system with arbitrary parameters. In this case we need the r.h.s of the system that will be solved: sng_system
	The vector of derivatives elements: jac
	The dimension of the system of equations: 9 
	A pointer to the arbitrary parameters of the system: (void*)&ode_param */
	gsl_odeiv2_system   ode_sys   = {sng_system, jac, 9, (void*)&ode_param};

/* GrowingMode is linear growing mode, interpolated on the grid. See cosmo.c for details (double GrowingMode(double z, double k))*/
	double D_in =  GrowingMode(1./amin-1.,1./Smoothing.Radius[ismooth]);

	double y[9] = {l1*D_in, l2*D_in, l3*D_in,
		          l1*D_in/(l1*D_in - 1.), l2*D_in/(l2*D_in - 1.), l3*D_in/(l3*D_in - 1.),
		          l1*D_in, l2*D_in, l3*D_in};

    /*---------------- Integration step --------------- */
	double olda = mya; 
	double oldlam = y[0]; 

	while (mya < amax) {
        /* gsl_odeiv2_evolve_apply advances the system from time t and position y using the stepping function selected before i.e. rkf45. 
		The new time and position are stored in t and y on output*/
        int status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &mya, amax, &hh, y);

        /* Correct integration check */
        if (status != GSL_SUCCESS) {
            printf("ERROR on task %d: integration of cosmological quantities failed\n", ThisTask);
            fflush(stdout);
            return -1;
        }

        /* ---------------- Update ellipsoid axis --------------- */
        if (y[0] >= 0.99999) {

            return olda + (1. - oldlam) * (mya - olda) / (y[0] - oldlam);

        }
    }

    /* In this case the ellipsoid does not collapse */
    return 0;
}

/*---------------------------------------- Calculation of b_c == growing mode at collapse time ---------------------------------------*/

double ell(const int ismooth,
	   const double l1,
	   const double l2,
	   const double l3)
{

#ifdef ELL_CLASSIC

    const double bc = ell_classic(ismooth, l1, l2, l3);

    if (bc > 0.0)
      return (1.0 + InverseGrowingMode(bc, ismooth));
    else
      return 0.0;

#endif // ELL_CLASSIC

#ifdef ELL_SNG

    const double bc = ell_sng(ismooth, l1, l2, l3);

    const double ret = ((bc > 0.0) ? (1.0 / bc) : 0.0);

    return ret;

#endif // ELL_SNG
}

/* ------------------------------------  Computation of collapse time i.e. F = 1 + z_collapse and variance ------------------------------------ */

/* Function: common_initialization */
/* Performed once by the host */
void common_initialization(const unsigned int size)
{
  for (unsigned int i=0 ; i<size ; i++)
    {
      /*----------- Common initialization ----------- */
      products[i].Rmax = -1;
      products[i].Fmax = (PRODFLOAT)-10.0;
            
      /*----------- Zel'dovich case ----------- */
      products[i].Vel[_x_] = (PRODFLOAT)0.0;
      products[i].Vel[_y_] = (PRODFLOAT)0.0;
      products[i].Vel[_z_] = (PRODFLOAT)0.0;

#ifdef TWO_LPT

      products[i].Vel_2LPT[_x_] = (PRODFLOAT)0.0;
      products[i].Vel_2LPT[_y_] = (PRODFLOAT)0.0;
      products[i].Vel_2LPT[_z_] = (PRODFLOAT)0.0;

#ifdef THREE_LPT

      products[i].Vel_3LPT_1[_x_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_1[_y_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_1[_z_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_2[_x_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_2[_y_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_2[_z_] = (PRODFLOAT)0.0;

#endif // THREE_LPT
#endif // TWO_LPT
    } // loop over MyGrids[0].total_local_size
  
  return;
}

int compute_collapse_times(int ismooth)
{
  /*----------------------------------------------------------------*/
  // Size of the products array
  const unsigned int total_size = MyGrids[0].total_local_size;
  
  if (!ismooth)
    {
      /*----------- Common initialization ----------- */
      common_initialization(total_size);
    }

  /* timing the main loop of 'compute_collapse_times' function */
  double cputmp = MPI_Wtime();
  /*-----------------------------------------------------------------------------------*/
  
  /*----------------- Calculation of variance, average, and collapse time -------------*/
  
  /* Local average and variance declaration */
  double local_average = 0.0, local_variance = 0.0;
  int all_fails = 0; 
  
  for (unsigned int index=0 ; index<total_size ; index++)
    {
      /* Computation of second derivatives of the potential i.e. the gravity Hessian */
      double diff_ten[6]; 
      for (int i=0 ; i<6 ; i++)
	{
	  diff_ten[i] = second_derivatives[0][i][index];
	}
        
      /* Computation of the variance of the linear density field */
      const double delta = (diff_ten[0] + diff_ten[1] + diff_ten[2]);
      local_average  += delta;
      local_variance += (delta * delta);

      /* Computation of the collapse time */
      double lambda1, lambda2, lambda3;
      int    fail;	
      /* inverse_collapse_time(funzione ell) */
      const double Fnew = inverse_collapse_time(ismooth, diff_ten, &lambda1, &lambda2, &lambda3, &fail);
      all_fails += fail;
      
      /* Updating collapse time */      
      const int update = (products[index].Fmax < Fnew);
      products[index].Rmax = (update ? ismooth : products[index].Rmax);
      products[index].Fmax = (update ? Fnew    : products[index].Fmax);

    } // target region  
   
  /* Fail check during computation of the inverse collapse time */
  /* If there were failures, an error message is printed and the function returns 1 */
  if (all_fails)
    {
      printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
      fflush(stdout);
      return 1;
    }
  
  /* Calculating obtained variance and avarage */		
  double global_variance = 0.0;
  double global_average  = 0.0;
    
  /* Calculates the global variance and average by reducing the values of local_variance and local_average across all MPI tasks */
  MPI_Reduce(&local_variance, &global_variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_average , &global_average , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* If the current task is the root task (task 0), it divides the global variance and average by the total number of data points (MyGrids[0].Ntotal) 
     to obtain the average values per data point */
  if (!ThisTask)
    {
      global_variance /= (double)MyGrids[0].Ntotal;
      global_average  /= (double)MyGrids[0].Ntotal;
    }
	
  /* Broadcasts the value of global_variance from the root task to all other tasks using MPI_Bcast */
  MPI_Bcast(&global_variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Stores the true variance*/
  Smoothing.TrueVariance[ismooth] = global_variance;

  /* CPU collapse time */	
  cputime.coll += (MPI_Wtime() - cputmp);
  
  return 0;

} // compute_collapse_times

/* ------------------------------------  Diagonalization of the potential Hessian ------------------------------------ */


double inverse_collapse_time(const int                     ismooth,
			     const double * const restrict dtensor,
				   double * const restrict x1,
				   double * const restrict x2,
				   double * const restrict x3,
				   int    * const restrict fail)
{  
  /* Local variables declaration */
  /* mu1, mu2 and mu3 are the principal invariants of the 3x3 tensor of second derivatives */
  /*  Performs diagonalization of the tensor by calculating the values of mu1, mu2, and mu3 */
  const double add[6] = {dtensor[0] * dtensor[0],
			 dtensor[1] * dtensor[1],
			 dtensor[2] * dtensor[2],
			 dtensor[3] * dtensor[3],
			 dtensor[4] * dtensor[4],
			 dtensor[5] * dtensor[5]};
  
  const double mu1   = (dtensor[0] + dtensor[1] + dtensor[2]);
  const double mu1_2 = (mu1 * mu1);
  const double mu2   = ((0.5 * mu1_2) - (0.5 * (add[0] + add[1] + add[2])) - (add[3] + add[4] + add[5]));
  
  /* mu3 calculation */
  const double mu3 = ((dtensor[0] * dtensor[1] * dtensor[2])       +
		      (2.0 * dtensor[3] * dtensor[4] * dtensor[5]) -
		      (dtensor[0] * add[5])                        -
		      (dtensor[1] * add[4])                        -
		      (dtensor[2] * add[3]));

  /* Check if the tensor is already diagonal */
  const double q = (mu1_2 - 3.0 * mu2) / 9.0;
  
  if (q == 0.0)
    {
      /* In this case the tensor is already diagonal */
      *x1 = dtensor[0];
      *x2 = dtensor[1];
      *x3 = dtensor[2];
    }
  else
    {  
      /* The solution has to be chosen as the smallest non-negative one between x1, x2 and x3 */
      const double r = -(((2.0 * mu1_2 * mu1) - (9.0 * mu1 * mu2) + (27.0 * mu3)) / 54.0);
        
      /* Fail check */
      if ((q * q * q) < (r * r) || (q < 0.0))
	{
	  *fail = 1;
	  return -10.0;
	}

      /* Calculating x1, x2, x3 solution in the same way as it done in ell_classic */
      const double sq = (2.0 * sqrt(q));
      const double t = acos(2.0 * r / q / sq);
      *x1 = ((-sq * cos(t * INV_3)) + (mu1 * INV_3));
      *x2 = ((-sq * cos((t + 2.0 * PI) * INV_3)) + (mu1 * INV_3));
      *x3 = ((-sq * cos((t + 4.0 * PI) * INV_3)) + (mu1 * INV_3));		
    }

  /* Ordering and inverse collapse time */
  ord(x1, x2, x3);
  
  /* Final computation of collapse time*/
  double t = ell(ismooth, *x1, *x2, *x3);
	

  *fail = 0;
  
  return t;
}

#define _MIN_(a,b) ((a) < (b) ? (a) : (b))
#define _MAX_(a,b) ((a) > (b) ? (a) : (b))

/* Orders a,b,c in decreasing order a>b>c */
void ord(double *const restrict a,
	 double *const restrict b,
	 double *const restrict c)
{
  const double hi = _MAX_(_MAX_(*a, *b), *c);
  const double lo = _MIN_(_MIN_(*a, *b), *c);
  *b = *a + *b + *c - lo - hi;
  *a = hi;
  *c = lo;

  return;
}
