/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2009

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file nestedsampler.cxx

*/

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

#include <boost/assign/list_of.hpp>

#include "nestedsampler.h"
#include "prior_sampler.h"
#include "nestedinitial.h"

namespace Minim {

#define INMIN_T(IN, L, M) if(IN->mon.trace_f && L <= IN->mon.trace_l) \
{ \
  IN->mon.trace_f(L, M, IN->mon.trace_d); };

  NestedS::NestedS(const inmin_nested_in * in,
		   const std::list<MCPoint> & start,
                   void *data,
		   unsigned seed):
    Zseq(1,0.0),
    Xseq(1,1.0),
    ps(new CSRMSSS(in, *this, g_ss())),
    initials(new InitialWorst()),
    in(in),
    cpoint(in->N),
    data(data),
    n_psample(100)
  {
    llPoint(in,
	    start,
	    ss,
            data);
    INMIN_T(in, 5, "Initialised the likelihoods of the starting set")
  }

  NestedS::~NestedS(void)
  {
  }

  void NestedS::reset(const std::list<MCPoint> &start,
		      CPriorSampler *cps )
  {
    ps.reset(cps);
    Zseq=boost::assign::list_of(0.0);
    Xseq=boost::assign::list_of(1.0);
    llPoint(in,
	    start,
	    ss,
            data);    
  }

  void NestedS::InitalS(NestedInitial *ins)
  {
    initials.reset(ins);
  }

  size_t NestedS::N(void) const
  {
    return ss.size();
  }

  double NestedS::sample(size_t j)
  {
    for (size_t i=0; i<j; ++i)
    {
      
      INMIN_T(in, 5, 
              (boost::format("Starting nested sampler iteration %i, have %i points in live set") % i % ss.size()).str().c_str());

      std::set<MCPoint>::iterator worst( --ss.end() );
      const double Llow=exp(-worst->ll);
      INMIN_T(in, 6, 
              (boost::format("Removed lowest point (lnLkl=%g)") % (-worst->ll)).str().c_str())

      const double X=exp(-((double)Xseq.size())/N());
      const double w=Xseq[Xseq.size()-1]-X;

      // Look for the next sample
      put(worst->p);
      const double newl = ps->advance((*initials)(*this).p,
                                      worst->ll,
				      n_psample);

      INMIN_T(in, 6, (boost::format("Adding new point to live set with lnLkl=%g") % newl).str().c_str());

      // Create new point
      MCPoint np;
      np.p.resize(in->N);
      get(np.p);
      np.ll=-newl;

      if (newl > 700)
      {
	//Reached maximum likelihood, exponential function and double
	//precision blow up after this!
	std::ostringstream msg;
	msg<<"Likelihood too high! Likelihood is:" << newl<<" at: "<<np;
	throw std::runtime_error(msg.str());
	//break;
      }

      // Is the new sample actually inside the contours of last?
      const bool better = -newl  < worst->ll;

      if (not better )
      {
        INMIN_T(in, 3, 
                "Could not find a better point so terminating early");
                
	// Note that this test is not definitive, since some
	// strategies will start from a point which not the worst
	// point and hence will return a "better" point event though
	// they haven not actually advanced their chain at all. See
	// below.
	break;
      }

      // Save the point about to be bumped off
      Zseq.push_back(Zseq[Zseq.size()-1] + Llow* w);
      Xseq.push_back(X);
      post.push_back(WPPoint(*worst, w));

      // Erase old point
      ss.erase(worst);
      
      std::pair<std::set<MCPoint>::iterator, bool> r=ss.insert(np);
      if(in->mon.accept_f)
        in->mon.accept_f(&np.p[0], np.ll, in->mon.accept_d);
      if (not r.second)
      {
        INMIN_T(in, 3, 
                "Could not insert a point because it has identical"
                "likelihood to an existing point. Can not contiue as we would have"
                "fewer points in the live set now.");
	break;
      }

    }
    return Zseq[Zseq.size()-1];
  }

  double NestedS::Z(void) const
  {
    return Zseq[Zseq.size()-1];
  }

  
  const std::list<WPPoint> & NestedS::g_post(void) const
  {
    return post;
  }
  
  const std::set<MCPoint> & NestedS::g_ss(void) const
  {
    return ss;
  }

  void llPoint(const inmin_nested_in * in,
	       const std::list<MCPoint> &lp,
	       std::set<MCPoint> &res,
               void *data)
  {
    for(std::list<MCPoint>::const_iterator i(lp.begin());
	i != lp.end();
	++i)
    {
      MCPoint p(i->p);
      // All likelihoods are on an inverted scale in these routines!
      p.ll=-1*in->fll(&p.p[0], data);
      res.insert(p);
    }
  }


  void printSS(const std::set<MCPoint> &ss)
  {
    for(std::set<MCPoint>::const_iterator i(ss.begin());
	i != ss.end();
	++i)
    {
      std::cout<<"p:";
      for(size_t j=0; j<i->p.size(); ++j)
	std::cout<<i->p[j]
		 <<",";
      std::cout<<i->ll<<",";
      std::cout<<std::endl;
    }
  }

}


