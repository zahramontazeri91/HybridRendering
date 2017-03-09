/*
   (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  
*/
/**  \file inmin_mon.h

   Monitoring of the minimisation/inference process
*/
#ifndef __INMIN_MON_H__
#define __INMIN_MON_H__

#include "stddef.h"

#ifdef __cplusplus
extern "C" {
#endif


  /** \brief Data to be passed to callbacks monitoring function
      evaluations
   */
  typedef struct  {

    /** Id of the execution that made this function evaluation, e.g.,
	proccess/thread/cpu/node ID. 
     */
    size_t execid;

    /** Parameters at which the function was evaluated. 
     */
    double *p;

  } inmin_mon_feval_d;


  /** \brief Type of pointer to callback for monitoring the function
      evaluations.

      \param d Information about current function evaluation

      \param state An opaque pointer to the state of the monitoring
      function
   */
  typedef  void (*inmin_mon_feval_f)(const inmin_mon_feval_d *d,
				     void *state);

  /** \brief Data to be passed to callbacks monitoring of the
      evaluation of Jacobians
   */
  typedef struct {
    /** Pointer to the jacobian matrix, size NxM
     */
    double *jac;
  } inmin_mon_jac_d;

  /** Type of pointer to callback for monitoring the evaluation sof
      the Jacobian matrix
      
      \param Structure containing information on the Jacobian that has
      been computed
      
      \param Opaque user supplied state
   */
  typedef  void (*inmin_mon_jac_f)(const inmin_mon_jac_d *d,
				     void *state);

  /** Callback for tracing algorithms
      
      \param level The tracin level to whic this message belongs

      \param msg The actual message
      
      \param state Pointer to an opaque state passed by the user
   */
  typedef  void (*inmin_mon_trace)(size_t level,
				   const char *msg,
				   void *state);


  

  /** \brief Monitoring definition to be passed to algorithms
   */
  typedef struct {

    /** Callback to monitor the function evaluations
     */
    inmin_mon_feval_f feval_f;
    /**  Opaque data to pass to the feval_f function
     */
    void * feval_d;

    /** Callback to monitor evaluation of Jacobians
     */
    inmin_mon_jac_f jac_f;
    /** Opaque data to pass to jac_f */
    void * jac_d;

    /** Callback for general tracing calls
     */
    inmin_mon_trace trace_f;
    /** Tracing level (0 is lowest)
     */
    size_t trace_l;
    /** Opaque data to pass trace_f
     */
    void * trace_d;

    /** Callback to monitor when a new point is accepted
     */
    void (*accept_f)(double *x, double ll, void *state) ;
    void *accept_d;

    /** Callback to monitor proposed points */
    void (*propose_p_f)(double *x, double ll, double p, void *state) ;
    void *propose_p_d;

  } inmin_mon_ctr;


  /** Simple ready made tracing function */
  void inmin_tracefn_stdout(size_t level,
                            const char *msg,
                            void *state);


#ifdef __cplusplus
}
#endif

#endif
