/*
   (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

   \file eigen_lm.cc

   Levenberg-Marquardt implementation using the Eigen library
*/   
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>
#include <boost/format.hpp>

#include "inmin_lm.h"
#include "eigen_utils.hpp"
#include "eigen_qo.hpp"

using namespace Eigen;

namespace inmin {

/// Class to carry out the evaluation of the supplied function
class FEval {

  inmin_f_res f;
  const size_t N;
  const size_t M;
  void *data;
  inmin_mon_ctr mon;

public:

  // Number of function evaluations
  size_t k;
  
  FEval(const inmin_lm_in *in,
	void *data):
    f(in->f),
    N(in->N),
    M(in->M),
    data(data),
    k(0)
  {
    memcpy((void *)&mon, 
	   (void *)&in->mon, 
	   sizeof(inmin_mon_ctr));
  }

  /** Evalue the function at x and store the result in res
   */
  void operator() (const VectorXd &x,
		   VectorXd &res)
  {
    f(&x[0], &res[0], data);
    ++k;

    if(mon.feval_f)
    {
      inmin_mon_feval_d md;
      md.execid=0;
      md.p=const_cast<double *>(&x[0]);
      mon.feval_f(&md, 
		  mon.feval_d);
    }
  }

};

/** Compute the Jacobian by a forward difference approximation
    
    \param x Compute the jacobian at this point
    \param fx is the value of functions at x
 */
void FwdDiffJac(const VectorXd &x,
		const VectorXd &fx,
		FEval &f,
		MatrixXd &B,
		double deltaa, 
		double deltar)
{
  VectorXd fxprime(fx.size());
  VectorXd xprime(x.size());
  for(size_t i=0; i<(size_t)x.size(); ++i)
  {
    const double delta=std::max(deltar*std::fabs(x[i]), deltaa);
    xprime=x;
    xprime[i]+=delta;
    f(xprime, fxprime);
    B.col(i)=(fxprime-fx)/delta;
  }
}

void LMCovar(const Eigen::MatrixXd &A,
	     Eigen::MatrixXd &covar,
	     double fnorm2,
	     size_t M)
{
  Eigen::MatrixXd Ap;
  size_t rank=PseudoInv(A, Ap);
  covar= Ap* (fnorm2/(M-rank));
}

}

#define INMIN_T(IN, L, M) if(IN->mon.trace_f && L <= IN->mon.trace_l) \
{ \
  IN->mon.trace_f(L, M, IN->mon.trace_d); };

using namespace inmin;

extern "C"  void inmin_lm_run_eigen(const inmin_lm_in *in,
				    const double *x0,
				    inmin_lm_out *out,
				    void *data)
{
  // Maximum mu before bailing out
  const double mumax =1e100;

  const inmin_lm_opts_s *opts;
  inmin_lm_opts_s defopts;
  if (in->opt.useopts)
  {
    opts=&(in->opt);
  }
  else
  {
    inmin_lm_opts_def(&defopts);
    opts=&defopts;
  }

  // Scaling for the damping factor
  double nu=2;

  // x is the current best guess
  // Defeat the supplied const to allow the creation of the view
  VectorXd x=Map<VectorXd>(const_cast<double *>(x0), 
			   in->N);
  // Store the result of f at x
  VectorXd fx(in->M);  

  // The approximate Jacobian 
  MatrixXd B(in->M, in->N);
  // Number of times we've approximated B by updates
  size_t Bapxs=0;

  FEval f(in, data);

  f(x, fx);

  // Compute the Jacobian directly
  FwdDiffJac(x, fx, f, B,
	     opts->dif_deltaa, 
	     opts->dif_deltar);

  
  MatrixXd A=B.transpose()*B;

  // The damping coefficient
  double mu=opts->tau*A.diagonal().maxCoeff();

  // The gradient
  VectorXd g=B.transpose()*fx;

  // Unless anything else is set, the stop is due to maximum number of
  // evaluations
  out->stop=LM_MAXFEV;

  // First check of gtol
  bool found=g.lpNorm<Infinity>()  < in->gtol;
  if (found)
  {
    out->stop=LM_GTOL;
    INMIN_T(in, 2,  "Stopping because maximum component of gradient is less than gtol");
  }
	   

  while( !found and f.k < in->maxfev)
  {
    if (mu  > mumax)
    {
      INMIN_T(in, 2,  "Stopping because damping is too large, can't make further progress");
      out->stop=LM_MUOVERFLOW;
      break;
    }
    // Damped A
    MatrixXd Ad=(A+mu*MatrixXd::Identity(A.rows(), A.cols()));
    VectorXd hlm;
    if (in->box)
    {  // Do constrained quadratic optimisation based on the damped A matrix
      VectorXd xnew;
      GQO_box(Ad, g-Ad*x, x, in->box, xnew);
      hlm=xnew-x;
    }
    else // Simple unconstrained solve
    {
      hlm=Ad.colPivHouseholderQr().solve(-1*g);
    }
      
    if(hlm.norm() < in->xtol*(x.norm()+in->xtol))
    {
      found=true;
      out->stop=LM_XTOL;
      INMIN_T(in, 2,  "Stopping because relative step size is smaller then xtol");
    }
    else
    {
      VectorXd xprime=x+hlm;
      VectorXd fxprime(fx.size());
      f(xprime, fxprime);
      
      const double Fimprov=(fx.squaredNorm()-fxprime.squaredNorm());
      const double LL=0.5*hlm.transpose()*(mu*hlm-g);
      const double sigma=Fimprov/LL;
      INMIN_T(in, 3, 
	      (boost::format("Trial step; improvement in sum of sq: %g, expected %g") 
	       % Fimprov % LL).str().c_str() );
      if (sigma>0 && Fimprov > 0)
      {
	INMIN_T(in, 2, (boost::format("Step accepted (mu: %g, sumsq: %g, sigma: %g)") 
			% mu % fxprime.squaredNorm() % sigma).str().c_str() );

	if (std::fabs(1 - fxprime.squaredNorm()/fx.squaredNorm()) < in->ftol)
	{
	  found=true;
	  out->stop=LM_FTOL;
	  INMIN_T(in, 2,  "Stopping because relative reduction in sum of squares is less then required");
	}
	if (Bapxs > in->N or nu > 16 )
	{
	  INMIN_T(in, 3, "Recomputing the entire Jacobian using forwad difference approxiation");
	  FwdDiffJac(xprime, 
		     fxprime, 
		     f, 
		     B,
		     opts->dif_deltaa, 
		     opts->dif_deltar);
	  Bapxs=0;
	}
	else
	{
	  INMIN_T(in, 3, "Updating the Jacobian using Broyden’s Rank One Update");
	  VectorXd u= (fxprime-fx-B*hlm)/hlm.squaredNorm();
	  B=B+u*hlm.transpose();
	  Bapxs++;
	}
	x=xprime;
	fx=fxprime;
	A=B.transpose()*B;
	g=B.transpose()*fx;
	mu=mu*std::max(opts->mu_scalemin, 1-std::pow(2*sigma-1, 3));
	nu=2;

	if (g.lpNorm<Infinity>()  < in->gtol)
	{
	  found=true;
	  out->stop=LM_GTOL;
	  INMIN_T(in, 2,  "Stopping because maximum component of gradient is less than gtol");
	}
      }
      else
      {
	INMIN_T(in, 2, (boost::format("Step rejected (mu: %g, sumsq: %g), increasing damping ") 
			% mu % fxprime.squaredNorm()).str().c_str() );
	mu=mu*nu;
	nu=2*nu;

	if (false)
	{
	  INMIN_T(in, 3, "Large increase in damping on updated Jacobian, doing full recomputation now ");
	  FwdDiffJac(x, 
		     fx, 
		     f, 
		     B,
		     opts->dif_deltaa, 
		     opts->dif_deltar);
	  Bapxs=0;
	}
      }
    }
    Map<VectorXd> xout(out->x, 
		       in->N);
    xout=x;

    if (out->covar)
    {
      Map<MatrixXd> covar(out->covar, 
			  in->N,
			  in->N);
      MatrixXd covartemp;
      LMCovar(A, covartemp, fx.squaredNorm(), in->M);
      covar=covartemp;
    }
  }
  
  
  
  
  
  
  
  
  
}
