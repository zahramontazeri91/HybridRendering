/**
   \file prior_sampler.hxx
   Bojan Nikolic <bojan@bnikolic.co.uk>, <b.nikolic@mrao.cam.ac.uk>

*/

#include <boost/bind.hpp>

#include "prior_sampler.h"

#include "mcpoint.h"
#include "markovchain.h"
#include "nestedsampler.h"

#define INMIN_T(L, M)  if(in->mon.trace_f && L <= in->mon.trace_l)     \
{ \
  in->mon.trace_f(L, M, in->mon.trace_d); };

namespace Minim
{
  
  CPriorSampler::CPriorSampler(const inmin_nested_in *in):
    in(in)
  {
  }

  CPriorSampler::~CPriorSampler()
  {
  }
  
  double likelihood(NestedS &s,
		    const std::vector<double> &x)
  {
    s.put(x);
    return s.llprob();
  }

  double prior(NestedS &s,
	       const std::vector<double> &x)
  {
    s.put(x);
    return s.pprob();
  }

  CSRMSSS::CSRMSSS(const inmin_nested_in *in,
                   NestedS &s,
		   const std::set<MCPoint> &ss):
    CPriorSampler(in),
    ss(ss),
    s(s)
  {
  }

  CSRMSSS::~CSRMSSS()
  {
  }

  bool rejImposPrior(const MCPoint2 &c, 
		     const MCPoint2 &p)
  {
    if (p.p-c.p > 800)
      {
      // exponential will evaluate to zero anyway in doble precision
      return false;
      }
    else
      return true;
  }

  void CSRMSSS::initChain(void)
  {
    MarkovChain::fx_t flkl=boost::bind(likelihood, 
				       boost::ref(s),
				       _1);
    
    MarkovChain::fx_t fprior=boost::bind(prior, 
					 boost::ref(s),
					 _1);
    std::vector<double> ic(ss.begin()->p.size());
    s.get(ic);

    c.reset(new ILklChain(ic,
			  flkl,
			  fprior,
			  constrPriorL,
			  rejImposPrior));

    nprop=0;
  }

  double CSRMSSS::advance(const std::vector<double> &ic,
                          double L,
			  size_t maxprop)
  {
    const double sf=0.1;

    if (not c) 
      initChain();

    const size_t n=c->n;


    //std::vector<double> ic(n);
    //s.get(ic);
    c->reset(ic,
	     L);

    std::vector<double> cv, eigvals, eigvects;
    omoment2(ss, cv);
    principalCV(cv, 
		eigvals, eigvects);

    for(size_t j=0; j<eigvals.size(); ++j)
      eigvals[j]= pow(eigvals[j],0.5)*sf;
      //eigvals[j]= eigvals[j]*sf;

    for(size_t i=0; i<maxprop; ++i)
    {
      std::vector<double> sigmas(n,0);
      sigmas[nprop%n]=eigvals[nprop%n];
      eigenProp(*c,
		//sigmas,
		eigvals,
		eigvects);
      if(in->mon.propose_p_f)
        in->mon.propose_p_f((double*)(& (c->gcx()[0])),
                            c->gcl(), 
                            -1, in->mon.propose_p_d);
      ++nprop;
    }
    s.put(c->gcx());

    // Note the convention is to return the negative value
    return -c->gcl();
  }

}

