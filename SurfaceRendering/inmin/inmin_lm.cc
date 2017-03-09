/*
   (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

*/

#include "inmin_lm.h"

extern "C"  void inmin_lm_run_pda(const inmin_lm_in *in,
				  const double *x,
				  inmin_lm_out *out,
				  void *data);

extern "C"  void inmin_lm_run_eigen(const inmin_lm_in *in,
				    const double *x,
				    inmin_lm_out *out,
				    void *data);

extern "C"  void inmin_lm_run(const inmin_lm_in *in,
			       const double *x,
			       inmin_lm_out *out,
			       void *data)
{
  // Switch to the Eigen implementation
  inmin_lm_run_eigen(in, x, out, data);
}

extern "C" void inmin_lm_opts_def(inmin_lm_opts_s *opt)
{
  opt->dif_deltaa=1e-6;
  opt->dif_deltar=1e-6;

  opt->tau=1e-3;
  opt->mu_scalemin=0.3333;
}



