/**
   \file prior_sampler.hxx
   Bojan Nikolic <bojan@bnikolic.co.uk>, <b.nikolic@mrao.cam.ac.uk>

   Sampler of the prior distribution 
*/

#ifndef _BNMIN1_PRIOR_SAMPLER_HXX__
#define _BNMIN1_PRIOR_SAMPLER_HXX__

#include <vector>
#include <set>

#include <boost/scoped_ptr.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "../inmin_nested.h"


#include "mcpoint.h"

namespace Minim
{
  // Forward declarations
  class MarkovChain;
  class InitPntChain;
  class ILklChain;
  class NestedS;

  /** \brief Constrained prior sampler
      
      Generate a new point at a position with probabilitiy
      proportional to prior distribution and sujbject to constraint
      that the likelihood is higher than supplied limit
   */
  class CPriorSampler
  {

  protected:

    const inmin_nested_in *in;

  public:
    
    // -------------- Construction/Destruction ---------------------

    /**
       \note We need separately the prior and the likelihood hence
       inheritance from indepenedentPriors
     */
    CPriorSampler(const inmin_nested_in *in);

    virtual ~CPriorSampler();

    // -------------- Public Interface -----------------------------

    /** \brief Advance the model to the now sample point
        
        \param ic Is the current point

	\param L The minimum likelihood of the new point

	\param maxprop Take no more than specified number of proposals
       
	\return the likelihood of the new current point
     */
    virtual double advance(const std::vector<double> &ic,
                           double L,
			   size_t maxprop) = 0;

  };

  /** Base the monte-carlo steps on the points in the live live set
     
   */
  class CSRMSSS:
    public CPriorSampler
  {

    boost::scoped_ptr<ILklChain> c;
    const std::set<MCPoint> &ss;

    void initChain(void);

    size_t nprop;

    NestedS &s;

  public:

    // -------------- Construction/Destruction ---------------------

    /**

     */
    CSRMSSS(const inmin_nested_in *in,
            NestedS &s,
	    const std::set<MCPoint> &ss);

    ~CSRMSSS();

    // -------------- Public Interface -----------------------------

    double advance(const std::vector<double> &ic,
                   double L,
		   size_t maxprop);

  };

  

}

#endif
