/*
   (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

*/

/** \file inmin_lm.h

   Levenberg-Marquardt least-squares fitting
*/   
#ifndef __INMIN_LM_H__
#define __INMIN_LM_H__

#include "stddef.h"

#include "inmin_mon.h"

#ifdef __cplusplus
extern "C" {
#endif

  /** Type of a function returning residuals.

      E.g., To be minimised in the least-squares sense
      
      \param x The parameters at which the function should be
      evaluated
      
      \param res Pointer to an array in which to store the computed
      residuals of the function
      
      \param data Opaque user-supplied pointer that is passed into the
      function

   */
  typedef  void (*inmin_f_res)(const double *x,
			       double *res,
			       void   *data);

  /** Options for parallel execution.
   */
  typedef struct
  {

    /** If set to 1, do parallel execution using the threads
	mechanism.
     */
    int doparallel;
    
    /** A function to copy the data passed to the function to be
	minimised. The new copy must be completely thread-independent
	(no shared data structures of any kind) from original.
     */
    void *  (* f_copy) (void *d);

    /** A function delete an extra copy of data created by f_copy
     */
    void  (* f_del) (void *d);

  } inmin_lm_parll_s;

  /** \brief Various tuning options for Levenberg-Marquardt
      minimisation

      Use the function inmin_lm_opts_def to set these to default
      values.
   */
  typedef struct
  {
    /** Overal switch for the options structure. If this is zero
    default options will be used, not the ones specified in this
    structure
    */
    int useopts;

    /** \brief Absolute delta to use when computing finite differences
     */
    double dif_deltaa;
    
    /** Relative delta to use when computing finite differences
     */
    double dif_deltar;

    /** Scaling for initial value of the damping factor */
    double tau;

    /** The maximum amuont to shrink the damping by */
    double mu_scalemin;

  } inmin_lm_opts_s;

  /** \brief Set the Levenberg-Marquardt options to default 
   */
  void inmin_lm_opts_def(inmin_lm_opts_s *opt);

  /** \brief Inputs to the Levenberg-Marquardt least-squares
      minimisation
   */
  typedef struct 
  {
    /** The number of free parameters */
    size_t N;
    /** Number of residuals to be minimised */
    size_t M;
    /**  The function to be minimised */
    inmin_f_res f;
    /** Target relative error in the sum of squares. Minimisation
	stops when the actual change in sum of squares between
	iterations is less than ftol. */
    double ftol;
    /** Target relative error in the values of the parameters
     */
    double xtol;
    /** Target orthogonality between the residual vector and the
	columns of the Jacobian */
    double gtol;
    /** Maximum number of function evaluations to use in the
	minimisation
     */
    size_t maxfev;
    /** \brief Array defining a box contraint for the optimisation
	
	NULL if no box constraint. Otherwise, it is an array of length
	N*2, where 2*i,2*i+1 elements are the lower and upper bounds
	of the i-th parameters.
    */
    double *box;

    /** Information about parallel minimisation. 
	
	\bug Parallelism is not currently implemented
    */
    inmin_lm_parll_s parll;

    /** Monitoring specification. Allows monitoring (to the standard
        output or otherwise) of various values and conditions in the
        algorithm.
    */
    inmin_mon_ctr mon;

    /** Options for tuning of the minimisation algorithm */
    inmin_lm_opts_s opt;

  } inmin_lm_in; 

  /** \brief Possible reasons for stop of the Levenberg-Marquardt
      algorithm
   */
  typedef enum {LM_FTOL, LM_GTOL, LM_XTOL, LM_MAXFEV, LM_MUOVERFLOW} inmin_lm_stop_e;

  /** \brief Output of the Levenberg-Marquardt algorithm
   */
  typedef struct 
  {
    /**  The parameters that minimise the residuals. This array must be
	 pre-allocated to size inmin_lm_in.N
    */
    double *x;

    /** The covariance matrix of the fit. This array must be
	pre-allocated to size inmin_lm_in.N * inmin_lm_in.N, or be
	NULL.
    */
    double *covar;

    /** Reason for termination
     */
    inmin_lm_stop_e stop;
  } inmin_lm_out; 

  /** Run the Levenberg-Marquardt least-squares minimisation.
      
     \param in Specification of minimisation problem and options

     \param x The initial position in parameter space to start the
     minimisation from

     \param out Pointer to structure in which to store the results of
     minimisation

     \param data Opaque pointer to pass to the function to be
     minimised
   */
  void inmin_lm_run(const inmin_lm_in *in,
		    const double *x,
		    inmin_lm_out *out,
		    void *data);



#ifdef __cplusplus
}
#endif

#endif
