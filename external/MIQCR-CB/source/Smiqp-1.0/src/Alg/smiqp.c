/****************************************************************************
 //  This file is part of the src of "SMIQP",
 //  Copyright (C) 2019  Amélie Lambert
 //  All Rights Reserved.
 // 
 //     CEDRIC - CNAM 
 //     292 rue saint martin
 //     F-75141 Paris Cedex 03
 //     France
 //      
 //     amelie.lambert@cnam.fr    http://cedric.cnam.fr/~lamberta 
 // 
 //  This code is published under the Eclipse Public License.
 // Author: Amélie Lambert
 //****************************************************************************/

#include<fcntl.h>
#include<stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include<stdlib.h>
#include<unistd.h>
#include<getopt.h>
#include<assert.h>
#include<string.h>
#include<fcntl.h>
#include<getopt.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <sys/timeb.h>
#include <string.h>

#include"utilities.h"
#include"quad_prog.h"
#include"in_out.h"
#include "liste.h"

/********************** Branch & Bound ****************************/
#ifdef HAVE_CPLEX
#include"miqcr_cplex.h"
#include"sbb_miqcr_cplex.h"
#include "sbb_boxqp_cplex.h"
#include "boxqp_cplex.h"
#endif

#ifdef HAVE_SCIP
#include"miqcr_scip.h"
#include"sbb_miqcr_scip.h"
#include"sbb_boxqp_scip.h"
#include "boxqp_scip.h"
#endif

#include"sbb_miqcr_ipopt.h"
#include"sbb_boxqp_ipopt.h"



/***************************  presolve *****************************/
#ifdef HAVE_CPLEX
#include "obbt.h"
#endif
/***************************  local sol *****************************/
#ifdef HAVE_CPLEX
#include"local_sol_cplex.h"
#endif

#ifdef HAVE_SCIP
#include"local_sol_scip.h"
#endif

#include"local_sol_ipopt.h"

/************************* Solver SDP *****************************/
#include"solver_sdp.h"
#include"solver_sdp_mixed.h"



/*************************************************************************************/
/**************************  binary, general integer, continuous, or mixed pb ********/
/*************************************************************************************/

/****************************************************************************************************/
/****************************** GLOBAL VARIABLES ****************************************************/
/****************************************************************************************************/

/* Variables for caraterizing the class of the instance to solve */
int type; /* 0 if binary, 1 if general integer, 2 if pure continuous, 3 if mixed */
int cont_lin=0; /* indicates if there is ony linear constraints */
int cont_quad=0; /* indicates if there is ony quadratic constraints */
int obj_quad=1;
int unconstrained=0;
int is_not_0_1=0 ;
int nb_var_cont; 
int nb_var_int;
int nb_var_01;
int nb_cont_01;
int nb_cont_sup_1;
int is_not_lb_0=0;
int is_bounded =1;


/*Variables for measuring the time */
long start_time;
long time_init_sol;
long time_obbt;
long time_sdp;
long time_bb;

double INFTY = 10e9;

/*variables for files */
 char *buf_in;
 char *buf_out;
 
/*Variables for the solution algorithms */

 Liste liste;
 double ub_bb;
 double lb_bb;
 double * best_x;
 long nodecount;
 double sol_sdp;

/*Structure for solving the semidefinite relaxation*/
 SDP psdp;

/*Sructures for solving the equivalent convex formulation or the convex relaxation*/
 MIQCP qp;
 Q_MIQCP qqp;
 C_MIQCP cqp;  

int IS_CPLEX=0;
int IS_SCIP=0;

extern int init_param();

extern double MAX_SOL_BB;
extern int NB_PASS_OBBT_INIT;
extern int OBBT;
extern int OBBT_QUAD;
extern int OBBT_CONT;

static void write_json_double_array(FILE *stream, const char *name,
                                    const double *values, int length,
                                    int trailing_comma)
{
  int i;
  fprintf(stream, "  \"%s\": [", name);
  for (i = 0; i < length; i++)
    {
      if (i > 0)
	fputs(", ", stream);
      fprintf(stream, "%.17g", values[i]);
    }
  fprintf(stream, "]%s\n", trailing_comma ? "," : "");
}

static void write_json_int_array(FILE *stream, const char *name,
                                 const int *values, int length,
                                 int trailing_comma)
{
  int i;
  fprintf(stream, "  \"%s\": [", name);
  for (i = 0; i < length; i++)
    {
      if (i > 0)
	fputs(", ", stream);
      fprintf(stream, "%d", values[i]);
    }
  fprintf(stream, "]%s\n", trailing_comma ? "," : "");
}

/* Export the pure-binary MIQCR model data after the authors' Phase 1.  The
   companion Python code reconstructs the same x/y/z/t formulation generated
   by miqcr_cplex.c, but solves it with Gurobi. */
static int export_binary_miqcr_phase2(const char *path)
{
  FILE *stream;

  if (qqp->nb_int != qqp->n)
    {
      fprintf(stderr,
              "Phase-2 export currently supports pure-binary problems only.\n");
      return 0;
    }

  stream = fopen(path, "w");
  if (stream == NULL)
    {
      perror("Could not open the Phase-2 export file");
      return 0;
    }

  fprintf(stream, "{\n");
  fprintf(stream, "  \"format\": \"smiqp-miqcr-binary-v1\",\n");
  fprintf(stream, "  \"n\": %d,\n", qqp->n);
  fprintf(stream, "  \"nb_int\": %d,\n", qqp->nb_int);
  fprintf(stream, "  \"m\": %d,\n", qqp->m);
  fprintf(stream, "  \"p\": %d,\n", qqp->p);
  fprintf(stream, "  \"mq\": %d,\n", qqp->mq);
  fprintf(stream, "  \"pq\": %d,\n", qqp->pq);
  fprintf(stream, "  \"nb_y\": %d,\n", qqp->nb_y);
  fprintf(stream, "  \"nb_z\": %d,\n", qqp->nb_z);
  fprintf(stream, "  \"nb_col\": %d,\n", qqp->nb_col);
  fprintf(stream, "  \"constant\": %.17g,\n", qqp->cons);
  fprintf(stream, "  \"new_lambda_min\": %.17g,\n",
          qqp->new_lambda_min);
  fprintf(stream, "  \"phase2_beta_tolerance\": %.17g,\n", EPS_BETA);
  write_json_double_array(stream, "q", qqp->q, qqp->n * qqp->n, 1);
  write_json_double_array(stream, "new_q", qqp->new_q,
                          qqp->n * qqp->n, 1);
  write_json_double_array(stream, "c", qqp->c, qqp->n, 1);
  write_json_double_array(stream, "lower_bounds", qqp->l, qqp->n, 1);
  write_json_double_array(stream, "upper_bounds", qqp->u, qqp->n, 1);
  write_json_double_array(stream, "a", qqp->a, qqp->m * qqp->n, 1);
  write_json_double_array(stream, "b", qqp->b, qqp->m, 1);
  write_json_double_array(stream, "d", qqp->d, qqp->p * qqp->n, 1);
  write_json_double_array(stream, "e", qqp->e, qqp->p, 1);
  write_json_double_array(stream, "aq", qqp->aq,
                          qqp->mq * (qqp->n + 1) * (qqp->n + 1), 1);
  write_json_double_array(stream, "bq", qqp->bq, qqp->mq, 1);
  write_json_double_array(stream, "alphaq", qqp->alphaq, qqp->mq, 1);
  write_json_double_array(stream, "dq", qqp->dq,
                          qqp->pq * (qqp->n + 1) * (qqp->n + 1), 1);
  write_json_double_array(stream, "eq", qqp->eq, qqp->pq, 1);
  write_json_double_array(stream, "alphabisq", qqp->alphabisq, qqp->pq, 1);
  write_json_double_array(stream, "beta", qqp->beta,
                          qqp->n * qqp->n, 1);
  write_json_int_array(stream, "y_mask", qqp->nb_var_y,
                       qqp->n * qqp->n, 0);
  fprintf(stream, "}\n");

  if (fclose(stream) != 0)
    {
      perror("Could not close the Phase-2 export file");
      return 0;
    }
  return 1;
}

/* Disk checkpoint runs in a child so model construction cannot alter the
   ongoing bundle state. Rename publishes only a fully written export. */
void checkpoint_binary_miqcr(void)
{
  const char *path=getenv("SMIQP_CHECKPOINT_PATH");
  static time_t last=0;
  time_t now=time(NULL);
  pid_t child;
  int status;
  if (path==NULL || (last!=0 && now-last<30)) return;
  fflush(NULL);
  child=fork();
  if (child==0)
    {
      int i,j;
      char *tmp=malloc(strlen(path)+5);
      if (tmp==NULL) _exit(1);
      sprintf(tmp,"%s.tmp",path);
      psdp->beta=alloc_matrix_d(psdp->n,psdp->n);
      for(i=0;i<psdp->n;i++) for(j=i;j<psdp->n;j++)
        {
          double v=(i==j) ? psdp->beta_diag[i] :
            psdp->beta_1[ij2k(i,j,psdp->n)]+psdp->beta_2[ij2k(i,j,psdp->n)]-
            psdp->beta_3[ij2k(i,j,psdp->n)]-psdp->beta_4[ij2k(i,j,psdp->n)];
          psdp->beta[ij2k(i,j,psdp->n)]=v;
          psdp->beta[ij2k(j,i,psdp->n)]=v;
        }
      EPS_BETA=0.0001; /* Same support tolerance as the experiment export. */
      create_q_miqcp();
      compute_new_q_miqcr();
      compute_new_lambda_min_miqcr();
      if (!export_binary_miqcr_phase2(tmp) || rename(tmp,path)!=0) _exit(1);
      _exit(0);
    }
  if (child<0) { fprintf(stderr,"CB checkpoint fork failed.\n"); return; }
  if (waitpid(child,&status,0)==child && WIFEXITED(status) && WEXITSTATUS(status)==0)
    { last=now; printf("CB disk checkpoint completed.\n"); fflush(stdout); }
  else fprintf(stderr,"CB checkpoint failed; prior checkpoint retained.\n");
}


int main (int argc, char **argv)
{
 start_time = time(NULL);

 
  int i;
  /*structure for spatial bb*/
  struct miqcp_bab mbab;
  int status;
  
  /*initialisation of best upper and lower bounds*/
  ub_bb=MAX_SOL_BB;
  lb_bb = -MAX_SOL_BB;

  /* check solvers availables*/
  #ifdef HAVE_CPLEX
  IS_CPLEX=1;
  #endif
  #ifdef HAVE_SCIP
  IS_SCIP=1;
  #endif
  
  /************************************************************************************************/
  /*************************  Reading instance file  **********************************************/
  /************************************************************************************************/
  int size_buf_in;
  int size_buf_out;
  int file_sol;
  if (argc >=2)
    {
      size_buf_in = strlen(argv[1]);
      buf_in = alloc_string(size_buf_in +1);
      strcpy(buf_in,argv[1]);
    }
 
  if (argc == 3)
    {
      size_buf_out = strlen(argv[2]);
      buf_out = alloc_string(size_buf_out +1);
      strcpy(buf_out,argv[2]);
      file_sol = open(buf_out, O_RDWR | O_CREAT | O_TRUNC ,S_IRWXU | S_IRWXG | S_IRWXO);
      status = dup2(file_sol,STDOUT_FILENO);
      close(file_sol);
      if (status ==-1)
	{
	  printf("\n Output file error \n");
	}
    
      }
  
  if (argc>3 || argc==1) 
   {
  	printf("\nError : you must give at least an instance file\nUsage\n./smiqp data.dat result.sol\nThe file result.sol is optional\n\n");
  	return 0;
    }
  
  /************************************************************************************************/
  /***************************** Load the problem *************************************************/
  /************************************************************************************************/
  status = read_file();
  if (status ==0)
    return 1;
  
  /************************************************************************************************/
  /***************************** Initialisation of parameters *************************************/
  /************************************************************************************************/

  status = init_param(); 

  if (status !=1)
    printf("\nError in parameter file\n");
  
  printf("############################################################################\n############################################################################\n\nStart solving.\n\n");
  printf("Reading file \"%s\"\n\n %d integer variables\n %d continuous variables\n %d variables\n %d linear equality constraints\n %d linear inequality constraints\n %d quadratic equality constraints\n %d quadratic inequality constraints\n\n",buf_in,qp->nb_int, qp->n - qp->nb_int, qp->n, qp->m,qp->p,qp->mq,qp->pq);
	  
  /************************************************************************************************/
  /********************** Determine the type of the problem to solve*******************************/
  /************************************************************************************************/
  
  /* type : 0 if binary, 1 if general integer, 2 if pure continuous, 3 if mixed */
  if (qp->nb_int == 0)
    type = 2;
  else
    {
      if (qp->nb_int != qp->n)
	type =3;
      else
	{
	  for(i=0;i<qp->n;i++)
	    if(qp->l[i]!=0 || qp->u[i]!=1)
	      {
		type=1;
		break;
	      }
	  if (i==qp->n)
	    type=0;
	}
	
    }

  /*determine the types of the constraints*/
  cont_lin=0;
  unconstrained=0;
  cont_quad=0;
  
  
  if ((qp->m>0) && (qp->mq - 1 == 0))
    {
      if (qp->pq == 0 && qp->p==0)
	cont_lin = 1;
      else
	if ((qp->p>0) && (qp->pq -(2*qp->p*qp->n) ==0))
	      cont_lin=1;
	else
	  if(qp->pq >0)
	    cont_quad=1;
    }
  else
    {
      if ((qp->p>0) && (qp->pq -(2*qp->p*qp->n) ==0) && ((qp->mq -1 == 0) || qp->mq==0))
	cont_lin=1;
      else
	if(qp->mq>0 || qp->pq >0)
	  cont_quad=1;
      
      if ((qp->m==0) && (qp->p==0) && (qp->mq  == 0) && (qp->pq ==0))
	unconstrained=1;
    }

    
  if   ((qp->p==0) && (qp->pq ==0))
    OBBT_CONT=0;
  is_not_0_1 = 0;
  for (i=0;i<qp->n;i++)
    if (qp->u[i] != 1 || qp->l[i]!=0)
      {
	is_not_0_1 =1;
	break;
      }

  is_not_lb_0 = 0;
  for (i=0;i<qp->n;i++)
   if (qp->l[i]!=0)
     {
       is_not_lb_0=1;
       break;
     }

   /************************************************************************************************/
  /************** Check if lower and upper bounds on variables are provided **************************/
  /************************************************************************************************/
  for(i=0;i<qp->n;i++)
    if(qp->u[i] >= INFTY || qp->l[i]<=-INFTY)
      is_bounded =0;

  if (IS_CPLEX ==0 && is_bounded ==0)
    {
      printf("\n At least one upper or lower bound on a variable is missing\nProvide lower and upper bounds on each variables or install Cplex\n");
      return 1;
    }
  
  /************************************************************************************************/
  /************** Compute upper and lower bounds for variables *************** **********************/
  /************************************************************************************************/

  #ifdef HAVE_CPLEX
  if (is_bounded == 0)
    {
      is_bounded=1;
      initialize_bounds_with_obbt();
      for(i=0;i<qp->n;i++)
	if(qp->u[i] >= INFTY || qp->l[i]<=-INFTY)
	  {
	    printf("\n At least one upper or lower bound on a variable is missing after OBBT : probably unbounded problem\n");
	    return 1;
	  }
    }
  #endif

  /************************************************************************************************/
  /************** Change bounds of x_i to 0 and 1 if pure continuous problem **********************/
  /************************************************************************************************/
  if (type==2 && is_not_0_1==1 && is_bounded==1)
    {
      printf("\nPresolve: Change bounds to 0 and 1\n");
       change_bounds();
    }

  /************************************************************************************************/
  /************** Change bounds of x_i to 0 and u_i - l_i if pure general integer problem **********************/
  /************************************************************************************************/
  if (type==1 && is_not_lb_0==1)
    {
      change_bounds_gen_int();
    }

  
 nb_var_cont=0; 
 nb_var_int=0;
 nb_var_01=0;
 nb_cont_01=0;
 nb_cont_sup_1=0;

  for(i=0;i<qp->nb_int;i++)
    if (qp->u[i] !=1)
      nb_var_int++;
    else
      nb_var_01++;

  for(i=qp->nb_int;i<qp->n;i++)
    {
      if (qp->u[i] ==1 && qp->l[i]==0)
	nb_cont_01++;
      if (qp->l[i]>=1)
	nb_cont_sup_1++;
    }

  nb_var_cont= qp->n - nb_var_int - nb_var_01;
  
  /************************************************************************************************/
  /****************************  Compute initial upper bound **************************************/
  /************************************************************************************************/
  int compute_local_sol = 0;
  best_x=alloc_vector_d(qp->n);
  if (type ==2 || type ==3)
    {
      printf("\nStarting computing initial local solution\n");
     
       if(cont_lin && type ==2)
	{
	  #ifdef HAVE_CPLEX
	  compute_local_sol_init();
	  if ( qp->local_sol != NULL)
	    {
	      ub_bb= qp->sol_adm ;
	      best_x =copy_vector_d(qp->local_sol,qp->n);
	      printf("\nCompute initial upper bound with cplex succeeded \n");
	      printf("\nLocal solution value (cplex) : %lf\n",ub_bb);
	      compute_local_sol=1;
	    }
	  else
	    printf("\nCompute initial upper bound with cplex failed \n");
          #endif
	  #ifdef HAVE_SCIP
	  if (compute_local_sol ==0)
	    {
	      if (qp->nb_int==0)
		status =compute_local_sol_ipopt_init();
	      else
		status = compute_local_sol_scip_init();
	      if ( status ==1)
		{
		  ub_bb= qp->sol_adm ;
		  best_x=copy_vector_d(qp->local_sol,qp->n);
		  if (qp->nb_int ==0)
		    {
		      printf("\nCompute initial upper bound with ipopt succeeded \n");
		      printf("\nLocal solution value (ipopt) : %lf\n",ub_bb);
		    }
		  else
		    {
		      printf("\nCompute initial upper bound with scip succeeded \n");
		      printf("\nLocal solution value (scip) : %lf\n",ub_bb);
		    }
		  compute_local_sol=1;
		 
		}
	      else
		printf("\nCompute initial upper bound with scip failed \n");		 
	    }
	  #endif
	  #ifndef HAVE_SCIP
	  if (compute_local_sol ==0 && qp->nb_int ==0)
	    {
	      status =compute_local_sol_ipopt_init();
	      if ( status ==1)
		{
		  ub_bb= qp->sol_adm ;
		  best_x=copy_vector_d(qp->local_sol,qp->n);
		  printf("\nCompute initial upper bound with ipopt succeeded \n");
		  printf("\nLocal solution value (ipopt) : %lf\n",ub_bb);
		  compute_local_sol=1;
		}
	      else
		printf("\nCompute initial upper bound with ipopt failed \n");	
	    }
	  #endif
	}
       else
	 {
         #ifdef HAVE_SCIP
	   if (qp->nb_int == 0)
	     status = compute_local_sol_ipopt_init();
	   else
	     status = compute_local_sol_scip_init();
	   if ( status ==1)
	     {
	       ub_bb= qp->sol_adm ;
	       best_x=copy_vector_d(qp->local_sol,qp->n);

	       if (qp->nb_int ==0)
		 {
		   printf("\nCompute initial upper bound with ipopt succeeded \n");
		   printf("\nLocal solution value (ipopt) : %lf\n",ub_bb);
		 }
	       else
		 {
		   printf("\nCompute initial upper bound with scip succeeded \n");
		   printf("\nLocal solution value (scip) : %lf\n",ub_bb);
		 }
	       compute_local_sol=1;
	       
	     }
	   else
	     printf("\nCompute initial upper bound with scip failed \n");		 
	   #endif
	   #ifndef HAVE_SCIP
	   if (compute_local_sol == 0  && qp->nb_int ==0)
	     {
	       status = compute_local_sol_ipopt_init();
	       if (status ==1)
		 {
		   ub_bb= qp->sol_adm ;
		   best_x=copy_vector_d(qp->local_sol,qp->n);
		   printf("\nCompute initial upper bound with ipopt succeeded \n");
		   printf("\nLocal solution value (ipopt) : %lf\n",ub_bb);
		      compute_local_sol=1;
		 }
	       else
		 printf("\nCompute initial upper bound with ipopt failed \n");	
	     }
	   #endif
	 }
       time_init_sol = time(NULL);
       if (status ==1)
	 printf("\nTime for computing initial upper bound: %ld \n", time_init_sol - start_time);
    }
  /*******************************************************************************/
  /************************** Reduce bounds with obbt ****************************/
  /*******************************************************************************/
  int proceed_obbt=0;
  
  if (unconstrained == 0 && type ==2 && IS_CPLEX == 1)
    proceed_obbt=1;

   #ifdef HAVE_CPLEX

  if(proceed_obbt==1 && OBBT==1)
    {
      if (OBBT_QUAD ==1)
	printf("\nStarting presolve with OBBT-quad \n");
      else
	printf("\nStarting presolve with OBBT\n");

      if (OBBT_CONT == 1)
	printf("\nContracting bounds on variables and inequality constraints\n");
      else
	printf("\nContracting bounds on variables\n");
      /*******************************************************************************/
      /********* Solve the sdp problem with conic bundle *****************************/
      /*******************************************************************************/
    
      
     
      create_sdp_mixed();
      if (OBBT_QUAD ==1)
	{
	  printf("Starting solving SDP-relaxation for OBBT-quad\n");
	  compute_alpha_beta_sbb_miqcr();
	  time_sdp = time(NULL);
	  printf("\nPre-processing time (solver sdp for obbt): %ld \n", time_sdp - time_init_sol);
	}
	  /*******************************************************************************/
	  /************************** Reduce bounds with obbt ****************************/
	  /*******************************************************************************/
      
	
      int nb_pass_var, nb_pass_cont;
      
      create_c_miqcp(2);
      if (OBBT_QUAD ==1)
	{ 
	  compute_new_q_sbb_miqcr();
	  compute_new_lambda_min_sbb_miqcr();
	}
      
      printf("\nStarting Presolve with OBBT and branch and Reduce : nb_pass_obbt_init max : %d\n",NB_PASS_OBBT_INIT);
      presolve_with_obbt(&nb_pass_var, &nb_pass_cont, ub_bb,NB_PASS_OBBT_INIT);
      printf("\n nb iterations for variable bounds reduction : %d\n iterations for constraint bounds reduction : %d\n",nb_pass_var,nb_pass_cont);
      time_obbt = time(NULL);
      printf("\nPre-processing time (OBBT): %.ld \n", time_obbt - time_init_sol);
	  
      
      /*************************************************************************************/
      /****************remove all non necessary constraints deduced from obbt***************/
      /*************************************************************************************/
      update_qp();
    }
  
  #endif

  /* Allow benchmarking the published MIQCR Phase 1 without requiring the
     optional CPLEX/SCIP Phase-2 interface. */
  if (getenv("SMIQP_PHASE1_ONLY") != NULL && type == 0 && unconstrained == 0)
    {
      const char *phase2_export;
      const char *phase2_beta_tolerance;

      create_sdp();
      compute_alpha_beta_miqcr();
      time_sdp = time(NULL);
      printf("Pre-processing time (sdp): %ld \n", time_sdp - start_time);

      phase2_export = getenv("SMIQP_EXPORT_PHASE2");
      if (phase2_export != NULL)
	{
	  phase2_beta_tolerance = getenv("SMIQP_PHASE2_BETA_TOLERANCE");
	  if (phase2_beta_tolerance != NULL)
	    {
	      char *endptr;
	      double requested_tolerance =
		strtod(phase2_beta_tolerance, &endptr);
	      if (endptr == phase2_beta_tolerance || *endptr != '\0' ||
		  requested_tolerance < 0)
		{
		  fprintf(stderr,
		          "Invalid SMIQP_PHASE2_BETA_TOLERANCE: %s\n",
		          phase2_beta_tolerance);
		  return 1;
		}
	      EPS_BETA = requested_tolerance;
	    }

	  create_q_miqcp();
	  compute_new_q_miqcr();
	  compute_new_lambda_min_miqcr();
	  if (!export_binary_miqcr_phase2(phase2_export))
	    return 1;
	  printf("Phase-2 beta support tolerance: %.12g\n", EPS_BETA);
	  printf("Phase-2 model data written to %s\n", phase2_export);
	}
      printf("Phase-1-only run complete.\n");
      return 0;
    }
  
  /*******************************************************************************/
  /********* Algorithms for binary  unconstrained quadratic programs *******/
  /*******************************************************************************/
    if(type == 0 && unconstrained == 1 && (IS_CPLEX == 1 || IS_SCIP ==1))
	{
	   create_sdp_mixed();
	  
	  compute_alpha_beta_sbb_miqcr();
	  time_sdp = time(NULL);

	  create_c_miqcp(3);

	  compute_new_q_boxqp();
	  compute_new_lambda_min_sbb_miqcr();

	  #ifdef HAVE_CPLEX
	  boxqp_cplex();
          #endif
	  
	  if (IS_CPLEX==0)
	    {
            #ifdef HAVE_SCIP
	    boxqp_scip();
	    #endif
	    }
	  
	  time_bb = time(NULL);
	  
	  printf("\nBest feasible solution value: %.5lf\n", ub_bb);
	  printf("\nBest feasible solution:\n");
	  for (i=0;i<qp->n;i++)
	    if ( !Zero_BB(best_x[i]))
	      printf("x[%d] = %.2lf\n", i, best_x[i]);
	  printf("\n");
	  printf("Pre-processing time (sdp): %ld \n", time_sdp - start_time);
	  printf("\nBranch-and-bound time: %ld \n", time_bb - time_sdp);
	  printf("Total time: %ld \n",time_bb - start_time);
	  printf("Total number of nodes: %ld \n\n", nodecount);
	  printf("Solving complete.\n\n############################################################################\n############################################################################\n\n");
	}
  




  /*******************************************************************************/
  /********* Algorithms for general integer quadratic programs *******************/
  /*******************************************************************************/
  // IS CPLEX ==1 type==1 IS_SCIP==1
    if( (type == 0 && unconstrained ==0 || type==1) && (IS_CPLEX == 1 || IS_SCIP ==1) )
    {
      create_sdp();
      
      compute_alpha_beta_miqcr();
      time_sdp = time(NULL);
            
      create_q_miqcp();
      compute_new_q_miqcr();
      compute_new_lambda_min_miqcr();
      #ifdef HAVE_CPLEX
	miqcr_cplex();
      #endif
	if(IS_CPLEX ==0)
	  #ifdef HAVE_SCIP
	  miqcr_scip();
	#endif
      time_bb = time(NULL);

      printf("\nBest feasible solution value: %.5lf\n", ub_bb);
      if (type==1 && is_not_lb_0==1)
	change_bounds_back_gen_int();
      printf("\nBest feasible solution:\n");
      for (i=0;i<qp->n;i++)
	if ( !Zero_BB(best_x[i]))
	  printf("x[%d] = %.2lf\n", i, best_x[i]);
      printf("\n");
      
      
      printf("Pre-processing time (sdp): %ld \n", time_sdp - start_time);
      printf("\nBranch-and-bound time: %ld \n", time_bb - time_sdp);
      printf("Total time: %ld \n", time_bb - start_time);
      printf("Total number of nodes: %ld \n\n", nodecount);
      printf("Solving complete.\n\n############################################################################\n############################################################################\n\n");
    }

  
  // type == 2 ou 3
  if(type == 2 || type == 3 || !(IS_CPLEX == 1 || IS_SCIP == 1))
    {
      
      //unconstrained == 1 et pure continu
      if (unconstrained ==1 && type ==2)
	{
	
	  create_sdp_mixed();
	  compute_alpha_beta_sbb_miqcr();
	  time_sdp = time(NULL);
	  
	  /* 1  for boxqp*/
	  create_c_miqcp(1);
	   
	  liste = init_liste();
	  
	  mbab=create_mbab();
	 
	  compute_new_q_sbb_miqcr();
	  compute_new_lambda_min_sbb_miqcr();
	  
	   #ifdef HAVE_CPLEX	  
	   sbb_boxqp_cplex(mbab);
	   #endif

	    if (IS_CPLEX==0){
	      #ifdef HAVE_SCIP
	      sbb_boxqp_scip(mbab);
	       #endif
	    }
	    
	    if (IS_SCIP ==0 && IS_CPLEX ==0)
	      sbb_boxqp_ipopt(mbab);
	  
	  time_bb = time(NULL);
	  printf("\nBest lowerbound: %.5lf\n", lb_bb);
	  printf("\nBest feasible solution value: %.5lf\n", ub_bb);
	  if (type==2 && is_not_0_1==1 && is_bounded==1)
	    change_bounds_back();
	  if (type==1 && is_not_lb_0==1)
	    change_bounds_back_gen_int();
	  printf("\nBest feasible solution:\n");
	  for (i=0;i<qp->n;i++)
	   if ( !Zero_BB(best_x[i]))
	      printf("x[%d] = %.2lf\n", i, best_x[i]);
	  printf("\n");

	  if (compute_local_sol ==1)
	    printf("\nTime for computing initial upper bound: %ld \n", time_init_sol - start_time);
	  else
	    time_init_sol=start_time;
	  
	  printf("Pre-processing time (sdp): %ld \n", time_sdp - time_init_sol);
	  printf("\nBranch-and-bound time: %ld \n", time_bb - time_sdp);
	  printf("Total time: %ld \n", time_bb - start_time);
	  printf("Total number of nodes: %ld \n\n", nodecount);
	  printf("Solving complete.\n\n############################################################################\n############################################################################\n\n");
	}
      else
	{
	  /*************************************************************************************/
	  /******************  solve  the sdp if necessary  *******************************/
	  /**************************************************************************************/
	  if (!(proceed_obbt ==1 && OBBT==1))
	    create_sdp_mixed();
	  if (OBBT_QUAD ==0 ||  (OBBT_QUAD ==1 && !(proceed_obbt ==1 && OBBT==1) ))
	    {
	
	      compute_alpha_beta_sbb_miqcr();
	      time_sdp = time(NULL);
	      if (getenv("SMIQP_PHASE1_ONLY") != NULL)
		{
		  printf("Pre-processing time (sdp): %ld \n", time_sdp - start_time);
		  printf("Phase-1-only run complete.\n");
		  return 0;
		}
	 
	    }
      
	  printf("\nStarting Branch-and-bound process\n");
	  
	  /* 2  for sbb_miqcr*/
	  create_c_miqcp(2);
	  liste = init_liste();
	  
	  mbab=create_mbab();
	  
	  nodecount=0;
	  compute_new_q_sbb_miqcr();
	 
	  compute_new_lambda_min_sbb_miqcr();
	
	  #ifdef HAVE_CPLEX	  
	    sbb_miqcr_cplex(mbab);
	  #endif
	    if (IS_CPLEX==0){
	    #ifdef HAVE_SCIP
	      sbb_miqcr_scip(mbab);
	      #endif
	    }
	    if (IS_SCIP==0 && IS_CPLEX==0)
	      sbb_miqcr_ipopt(mbab);
	  
	  
	  
	  time_bb = time(NULL);
      
	  printf("\nBest lowerbound: %.5lf\n", lb_bb);
	  printf("\nBest feasible solution value: %.5lf\n", ub_bb);
	  if (type==2 && is_not_0_1==1 && is_bounded==1)
	    change_bounds_back();
	  if (type==1 && is_not_lb_0==1)
	    change_bounds_back_gen_int();
	  printf("\nBest feasible solution:\n");
	  for (i=0;i<qp->n;i++)
	   if ( !Zero_BB(best_x[i]))
	      printf("x[%d] = %.2lf\n", i, best_x[i]);
	  printf("\n");

	  if (compute_local_sol ==1)
	    printf("\nTime for computing initial upper bound: %ld \n", time_init_sol - start_time);
	  else
	    time_init_sol=start_time;
	  
	  if (proceed_obbt ==1 && OBBT==1 && OBBT_QUAD==1)
	    {
	      printf("\nPre-processing time (solver sdp for obbt): %ld \n", time_obbt - time_init_sol); 
	      printf("Pre-processing time (sdp): %ld \n", time_sdp - time_obbt);
	    }
	  else
	    printf("Pre-processing time (sdp): %ld \n", time_sdp - time_init_sol); 
	  printf("\nBranch-and-bound time: %ld \n", time_bb - time_sdp);
	  printf("Total time: %ld \n", time_bb - start_time);
	  printf("Total number of nodes: %ld \n\n", nodecount);
	  printf("Solving complete.\n\n############################################################################\n############################################################################\n\n");
	}
    }

  return 1;
}
