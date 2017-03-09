/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2009

   LICENSE: GNU General Public License V2  

   \file mcpoint.cc
*/

#include <cmath>
#include <algorithm>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <boost/foreach.hpp>

#include "mcpoint.h"

namespace Minim {
  
  MCPoint::MCPoint(void):
    p(0),
    ll(-9999),
    fval(0)
  {
  }
  
  MCPoint::MCPoint(const std::vector<double> &p):
    p(p),
    ll(-9999),
    fval(0)
  {
  }  

  MCPoint::MCPoint(size_t np):
    p(np),
    ll(-9999),
    fval(0)
  {
  }

  MCPoint::MCPoint(const MCPoint &other):
    p(other.p),
    ll(other.ll),
    fval(other.fval)
  {
  }

  MCPoint &MCPoint::operator=(const MCPoint &other)
  {
    p=other.p;
    ll=other.ll;
    fval=other.fval;
    return *this;
  }

  std::ostream & operator<<(std::ostream &o, const MCPoint &p)
  {
    o<<"Likelihood: "<<p.ll
     <<" at: ";
    BOOST_FOREACH(double x, p.p)
      o<<x<<",";
    return o;
  }

  void moment1(const std::list<WPPoint> &l,
	       std::vector<double> &res)
  {
    const size_t n=l.begin()->p.size();
    res=std::vector<double>(n, 0.0);
    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	  res[j]+= (i->p[j] * i->w * exp(- i->ll));
      }
    }
  }

  void moment1(const std::list<WPPoint> &l,
	       double Z,
	       std::vector<double> &res)
  {
    moment1(l,res);
    for(size_t j=0; j<res.size(); ++j)
      res[j] /= Z;
  }

  void moment2(const std::list<WPPoint> &l,
	       const std::vector<double> &m1,
	       std::vector<double> &res)
  {
    const size_t n=m1.size();
    res=std::vector<double>(n, 0.0);
    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	res[j]+= ( pow(i->p[j]-m1[j],2.0) * i->w * exp(- i->ll));
      }
    }
  }

  void moment2(const std::list<WPPoint> &l,
	       const std::vector<double> &m1,
	       double Z,
	       std::vector<double> &res)
  {
    moment2(l, m1, res);
    for(size_t j=0; j<res.size(); ++j)
      res[j] /= Z;    
  }

  void moment1(const std::set<MCPoint> &s,
	       std::vector<double> &res)
  {
    const size_t n=s.begin()->p.size();
    res=std::vector<double>(n, 0.0);

    size_t N=0;
    for(std::set<MCPoint>::const_iterator i=s.begin();
	i!= s.end();
	++i)
    {
      if(i->p.size() != n)
      {
        abort();
      }
      for (size_t j=0; j<n; ++j)
      {
	res[j]+= (i->p[j]);
      }
      ++N;
    }
    
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]/=N;
    }
  }

  void moment2(const std::set<MCPoint> &s,
	       const std::vector<double> &m1,
	       std::vector<double> &res)
  {
    const size_t n=m1.size();
    res=std::vector<double>(n, 0.0);

    size_t N=0;
    for(std::set<MCPoint>::const_iterator i=s.begin();
	i!= s.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	res[j]+= pow(i->p[j]-m1[j], 2);
      }
      ++N;
    }
    
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]/=N;
    }
  }

  void omoment2(const std::set<MCPoint> &s,
		const std::vector<double> &m1,
		std::vector<double> &res)
  {
    const size_t n=m1.size();
    res=std::vector<double>(n*n, 0.0);

    size_t N=0;
    for(std::set<MCPoint>::const_iterator i=s.begin();
	i!= s.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	for(size_t k=0; k<n; ++k)
	{
	  res[j*n+k] += (i->p[j]-m1[j])*(i->p[k]-m1[k]);
	}
      }
      ++N;
    }
    
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]/=N;
    }

  }
  void omoment2(const std::set<MCPoint> &s,
		std::vector<double> &res)
  {
    std::vector<double> m1;
    moment1(s, m1);
    omoment2(s, m1, res);
  }

  void StdDev(const std::set<MCPoint> &s,
	      std::vector<double> &res)
  {
    std::vector<double> m1, m2;
    moment1(s, m1);
    moment2(s, m1, m2);
    res.resize(m2.size());
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]=pow(m2[j],0.5);
    }
  }

  void principalCV(const std::vector<double> &cv,
		   std::vector<double> &eigvals,
		   std::vector<double> &eigvects)
  {
   using namespace Eigen;
   const size_t n=sqrt(cv.size());
   Map<MatrixXd> CV(const_cast<double *>(&cv[0]), n, n);

   SelfAdjointEigenSolver<MatrixXd> eig(CV);

   eigvals.resize(n);
   eigvects.resize(n*n);
   for(size_t j=0; j<n; ++j)
   {
      eigvals[j]=eig.eigenvalues()[j];
      for(size_t i=0; i<n; ++i)
      {
        eigvects[j*n+i]= eig.eigenvectors()(i,j);
      }
    }

   //std::cout<<eig.eigenvalues()<<std::endl;
   //std::cout<<eig.eigenvectors()<<std::endl;

  }

  void postHist(const std::list<WPPoint> &l,
		double Z,
		const std::vector<double> &low,
		const std::vector<double> &high,
		size_t nbins,
		std::vector<double> &res)
  {
    const size_t ndim=low.size();

    res.resize(pow(nbins,ndim));
    std::fill(res.begin(), res.end(), 0.0);
    

    std::vector<double> deltas(ndim);
    for(size_t i=0; i<ndim; ++i)
    {
      deltas[i]=(high[i]-low[i])/nbins;
    }

    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      bool inside=true;
      size_t k=0;
      for (size_t j=0; j<ndim; ++j)
      {
	int dimi = int((i->p[j]-low[j])/deltas[j]);
	if (dimi >= 0 and dimi < (int)nbins)
	{
	  k+= dimi * pow(nbins, ndim-j-1);
	}
	else
	{
	  inside=false;
	}
      }
      if (inside)
      {
	res[k]+= i->w * exp(- i->ll);
      }
    }
  }

  void marginHist(const std::list<WPPoint> &l,
		  size_t pi,
		  double Z,
		  double low,
		  double high,
		  size_t nbins,
		  std::vector<double> &res)
  {
    res.resize(nbins);
    std::fill(res.begin(), res.end(), 
	      0.0);

    const double d=(high-low)/nbins;
    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      int k=int((i->p[pi]-low)/d);
      if (k > 0 and k < (int)nbins)
      {
	res[k]+= i->w * exp(- i->ll);
      }
    }

    for(size_t i=0; i<res.size(); ++i)
    {
      res[i]/=Z;
    }
  }

  void marginHist2D(const std::list<WPPoint> &l,
		    double Z,
		    size_t i,
		    double ilow,
		    double ihigh,
		    size_t j,
		    double jlow,
		    double jhigh,
		    size_t nbins,
		    std::vector<double> &res)
  {
    // Two dimensions only
    res.resize(pow(nbins,2));
    std::fill(res.begin(), res.end(), 
	      0.0);
    const double idelta=(ihigh-ilow)/nbins;
    const double jdelta=(jhigh-jlow)/nbins;
    
    for(std::list<WPPoint>::const_iterator p=l.begin();
	p!= l.end();
	++p)
    {
      
      int dimi = int((p->p[i]-ilow)/idelta);
      int dimj = int((p->p[j]-jlow)/jdelta);
      
      if (dimi >= 0 and   dimi<((int)nbins)  and   dimj >= 0 and  dimj < ((int)nbins))
      {
	const size_t k= dimi*nbins + dimj;
	res[k]+= p->w * exp(- p->ll);
      }
      
    }
  }
}


