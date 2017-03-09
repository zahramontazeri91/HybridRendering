/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
   (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

*/

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

#include <Eigen/Core>

#include "inmin_nested.h"

#include "bnmin1/nestedsampler.h"
#include "bnmin1/mcpoint.h"
#include "bnmin1/nestedinitial.h"



extern "C" void inmin_nested_run(size_t n,
                                 const double *ss,
                                 size_t nss,
                                 const inmin_nested_in *in,
                                 inmin_nested_out *out,
                                 void *data)
{
  using namespace Minim;
  std::list<MCPoint>  start;
  for(size_t i=0; i<nss; ++i)
    {
      MCPoint p(in->N);
      for(size_t j=0; j<in->N; ++j)
        p.p[j]=ss[i*(in->N) + j];
      start.push_back(p);
    }
  
  NestedS sampler(in,
                  start,
                  data,
                  43);
  sampler.InitalS(new Minim::InitialRandom(nss));
  //sampler.InitalS(new InitialWorst());
  sampler.sample(n);

  out->N=in->N;
  out->ns=sampler.g_post().size();
  out->xv=(double *)malloc(sizeof(double)*out->ns * in->N);
  out->llv=(double *)malloc(sizeof(double)*out->ns);
  out->wv=(double *)malloc(sizeof(double)*out->ns);

  std::list<WPPoint>::const_iterator ii=sampler.g_post().begin();
  for(size_t i=0; i<out->ns; ++i)
      {
          for(size_t j=0; j<in->N; ++j)
              out->xv[i*in->N+j]=ii->p[j];
          out->llv[i]=-1*ii->ll;
          out->wv[i]=ii->w;
          ++ii;
      }
  
  
}

extern "C"   void inmin_flatbox_ss(const double *box,
                                     size_t N,
                                     double *res,
                                     size_t NS,
                                     unsigned seed)
{

   boost::mt19937 rng(seed);
   boost::uniform_01<boost::mt19937> zo(rng);

   for(size_t i=0; i<NS; ++i)
   {
       for(size_t j=0; j<N; ++j)
           res[i*N+j]= box[2*j] + (box[2*j+1]-box[2*j])*zo();
   }
}

extern "C"    double inmin_evidence_rect(size_t n,
                                         const double *llv,
                                         const double *wv)
{
    using namespace Eigen;
    VectorXd ll=Map<VectorXd>(const_cast<double *>(llv), 
                              n);
    VectorXd w=Map<VectorXd>(const_cast<double *>(wv), 
                             n);
    return (ll.array().exp()* w.array()).sum();
}

extern "C"  double inmin_boxprior_f(const double *x,
                                    inmin_boxprior_data *data)
{
  for(size_t i=0; i<data->N; ++i)
    if ( x[i] < data->box[2*i] || x[i] > data->box[2*i+1])
        return -1e9;
  return 0;
}

extern "C" void inmin_fan_histogram(const inmin_nested_out *nested_res,
                                    void (*yfn)(const double *x,
                                                size_t nbins,
                                                double *y,
                                                const double *p,
                                                void   *data),
                                    void *data, 
                                    double xlow, double xhigh,
                                    double ylow, double yhigh,
                                    size_t nbins,
                                    double *res)
{
    const double xdelt= (xhigh-xlow)/(nbins+1);
    const double ydelt= (yhigh-ylow)/(nbins+1);
    std::vector<double> y(nbins);
    std::vector<double> x(nbins);
    for(size_t i=0; i<x.size(); ++i)
        {
            x[i]=xlow+(i+0.5)*xdelt;
            y[i]=0;
        }
    memset(res, 0, sizeof(*res)*nbins*nbins);
    for(size_t i=0; i<nested_res->ns; ++i)
    {
        yfn(&x[0], nbins, &y[0],
            &nested_res->xv[nested_res->N*i],
            data);
        const double dl=exp(nested_res->llv[i])*nested_res->wv[i];
        for(size_t j=0; j<y.size(); ++j)
            {
                int yi= (y[j]-ylow)/ydelt +0.5;
                if (yi<0 || yi>= nbins)
                    continue;
                res[yi*nbins+j] +=dl;
            }
    }
}
