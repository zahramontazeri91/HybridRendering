/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
   (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

*/

/** \file inmin_nested.h

   Nested sampler Bayesian analysis
*/   
#ifndef __INMIN_NESTED_H__
#define __INMIN_NESTED_H__

#include "stddef.h"
#include "inmin_mon.h"

#ifdef __cplusplus
extern "C" {
#endif

    /** Likelihood function type declaration
        
        \param x The parameters at which the likelihood should be
        evaluated
        
        \param data Opaque user-supplied pointer that is passed into the
        function
        
        \return The log-likelihood of the point x
    */
    typedef  double (*inmin_f_ll)(const double *x,
                                  void   *data);

   /** \brief Inputs to the Nested sampling algorithm    */
   typedef struct 
   {
    /** The number of free parameters */
    size_t N;

    /**  The likelihood function. The returned value should be
         log-likelihood of the point passed as first parameter to
         the function */
    inmin_f_ll fll;

    /**  The prior function. The returned value is the log-prior
         probability of the supplied point.*/
    inmin_f_ll fprior;

    /** Monitoring specification. Allows monitoring (to the standard
        output or otherwise) of various values and conditions in the
        algorithm.
    */
    inmin_mon_ctr mon;

   } inmin_nested_in; 

    /** \brief Output of the Nested sampling algorithm
     */
    typedef struct 
    {
        /** Number of parameters, repeated here to make this
            structure complete */
        size_t N;

        /** Number of samples that were made during the run
         */
        size_t ns;

        /** Array containing the positions of sample points. Memory is
            allocated by the nested sampler but must be free()d by
            the user. Size is N*ns.
        */
        double *xv;

        /** Array containing the likelihood values. Memory is
            allocated by the nested sampler but must be free()d by the
            user */
        double *llv;

        /** Array containing the weights. Memory is allocated by the
            nested sampler but must be free()d by the user */
        double *wv;
    } inmin_nested_out; 
    
    /** Advance the nested sampler

        \param n Maximum number of samples to make

        \param ss Starting set (total size of array must be (in->N)*nss)

        \param nss Number of points in the starting set
        
        \param in Input data structure

        \param out Output data structure
        
        \param data Opaque data to pass to the likelihood and prior
        functions 

     */
    void inmin_nested_run(size_t n,
                          const double *ss,
                          size_t nss,
                          const inmin_nested_in *in,
                          inmin_nested_out *out,
                          void *data);

    /** Generate a starting set for a flat box prior

        \param box Array of lenngth 2*N defining the prior box by
        pairs of values bounding each parameter, i.e., [param0_low,
        param0_high, param1_low, param1_high,....]
        
        \param N number of parameters

        \param res Output array (size NS*N) into which the starting
        set is stored

        \param NS number of points in the starting set
        
        \param seed Seed for the random number generator
     */
    void inmin_flatbox_ss(const double *box,
                          size_t N,
                          double *res,
                          size_t NS,
                          unsigned seed);

    /** Evaluate the evidence using simple rectangular integration

        \param n Number of points in the sequence of likelihood/weight
        values
        
        \param llv Array (length n) of likelihood values (e.g.,
        inmin_nested_out.llv)

        \param wv Array (length n) of weight values (e.g.,
        inmin_nested_out.wv)
     */
    double inmin_evidence_rect(size_t n,
                               const double *llv,
                               const double *wv);
    
    typedef struct
    {
        /** Number of parameters */
        size_t N;
        /** Array of size 2*N containing the box */
        double *box;
            
    } inmin_boxprior_data;

    /** A function which can be used as a prior function for box-like
        priors
        
        \param data The opaque data should be a inmin_boxprior_data
        structure. Note that currently nested sampling interface
        passes the same data to both likelihood and prior function so
        another layer of inderiction is needed if likelihood function
        needs its own data too.
     */
    double inmin_boxprior_f(const double *x,
                            inmin_boxprior_data *data);

    /** Compute a fan diagram histogram of a one-dimensional function
        model

        \param yfn Pointer to the model function, which takes
        following values:
           const double *x: Array of abcisse values to evaluate the
           function at
           size_t nbins: Number of abcissa values (bins)
           double *y: Array to store the function value in (pre-allocated)
           const double *p: Parameter set to evalute the function at
           void   *data: Opaque data to pass to the function

        \param data Opaque data to pass the the likelihood function

        \param xlow, xhigh, ylow, yhigh (range of the axes, from
        bottom of lowest bin to top of highest bin)

        \param nbins Number of bins to use (same in both x and y
        directions)

        \param res Output the result here. Must be a nbins X nbins
        pre-allocated array
     */ 
    void inmin_fan_histogram(const inmin_nested_out *nested_res,
                             void (*yfn)(const double *x,
                                         size_t nbins,
                                         double *y,
                                         const double *p,
                                         void   *data),
                             void *data, 
                             double xlow, double xhigh,
                             double ylow, double yhigh,
                             size_t nbins,
                             double *res);


#ifdef __cplusplus
}
#endif


#endif
