/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
/*
   (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

*/

/** \file inmin_compile.h

   Compile functions for minimisation
*/
#ifndef __INMIN_COMPILE_H__

#include "inmin_lm.h"

#ifdef __cplusplus
extern "C" {
#endif

    /** Functions compiled by inmin_s1d_fn_c expect this structure to
        be passed to their opaque data pointer.
     */
    typedef struct {
        /// Number of residual data points
        unsigned M;
        /// Pointer to array of independent variable values
        const double *xa;
        /// Pointer to array of dependent variable values (e.g., observations)
        const double *ya;
    } inmin_s1d_data;

  /** \brief Compile a simple one dimensional fitting function written
      in C-like syntax
   */
  inmin_f_res inmin_s1d_fn_c(const char *f);
  

#ifdef __cplusplus
}
#endif

#endif
