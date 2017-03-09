/*
   (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

   \file eigen_qo.cc

   Quadratic optimisation using the Eigen library
*/   

#include <stdexcept>
#include <algorithm>
#include <map>
#include <iostream>
#include <limits>

#include <Eigen/QR>
#include "eigen_qo.hpp"
#include "eigen_utils.hpp"

using namespace Eigen;

namespace inmin {

   int BQO_raw(const MatrixXd &H,
               const VectorXd &g,
               const MatrixXd &A,
               const VectorXd &b,
               VectorXd &xout,
               VectorXd &lambda)
   {

       Eigen::MatrixXd HInv;
       size_t HInvRank=PseudoInv(H, HInv);

       if (A.cols() != 0)
         {
           lambda=(A.transpose()*HInv*A).colPivHouseholderQr().solve(b+A.transpose()*HInv*g);
           xout=HInv*(A*lambda - g);
         }
       else
         {
           lambda*=0;
           xout=HInv*(-g);
         }
       return 0;
   }

   std::set<size_t> Box_broken(const Eigen::VectorXd &x,
                               const double *box)
   {
     std::set<size_t> res;
     const double epsf=2*std::numeric_limits<double>::epsilon();
     for(size_t i=0; i<x.size(); ++i)
       {
         if (x[i]< box[2*i] - (std::abs(epsf*box[2*i]) +epsf) )
           res.insert(2*i);
         if (x[i]> box[2*i+1] + (std::abs(epsf*box[2*i+1])+epsf))
           res.insert(2*i+1);
       }
     return res;
   }

  std::set<size_t> Box_BestFeasible(const Eigen::VectorXd &x0,
                                    const Eigen::VectorXd &xeq,
                                    const double *box,
                                    const std::set<size_t> &newConstr,
                                    Eigen::VectorXd &xout)
  {
    std::set<size_t> res;
    xout=xeq;

    // Map from t-value to index of constraint
    std::map<double, size_t> tmap;
    
    for(std::set<size_t>::const_iterator i=newConstr.begin();
        i != newConstr.end();
        ++i)
      {
        size_t ci=*i;
        const double t=(box[ci]-x0[ci/2])/(xeq[ci/2]-x0[ci/2]);
        tmap.insert(std::make_pair(t, ci));
      }
    
    const double tmin=tmap.begin()->first;
    const double ci=tmap.begin()->second;
    res.insert(ci);
    xout=x0+ (xeq-x0)* tmin;
    return res;
  }

  void Box_MakeConstrA(const std::set<size_t> &active,
                       const double *box,
                       size_t n,
                       Eigen::MatrixXd &A,
                       Eigen::VectorXd &b)
  {
    A=MatrixXd(n, active.size());
    b=VectorXd(active.size());
    A*=0;
    b*=0;
    size_t j=0;
    for(std::set<size_t>::const_iterator i=active.begin();
        i != active.end();
        ++i)
      {
        size_t ci=*i;
        A(ci/2, j)=1.0;
        b(j)=box[ci];
        ++j;
      }
  }

  std::set<size_t> negativeIndices(const Eigen::VectorXd &x)
  {
    std::set<size_t> res;
    for(size_t i=0; i<x.size(); ++i)
      if(x[i]<0)
        res.insert(i);
    return res;

  }

  int GQO_box(const Eigen::MatrixXd &H,
              const Eigen::VectorXd &g,
              const Eigen::VectorXd &x0,
              const double *box,
              Eigen::VectorXd &xout)
  {
    std::set<size_t> active;
    Eigen::VectorXd x=x0;
    while(1)
      {
          Eigen::MatrixXd A;
          Eigen::VectorXd b;
          Box_MakeConstrA(active, box, g.size(), A, b);
          // Solve the "current" basic quadratic optimisation problem
          Eigen::VectorXd xeq, lambda;
          BQO_raw(H, g, A, b, xeq, lambda);

          // Reverse the sign of lambda values which correspond to upper
          // bounds
          size_t j=0;
          for(std::set<size_t>::const_iterator i=active.begin();
              i != active.end();
              ++i)
            {
              size_t ci=*i;
              if (ci % 2 ) // Upper bound, need to reverse lambda
                lambda[j]*=-1;
              ++j;
            }

          std::set<size_t> v=Box_broken(xeq, box);
          std::set<size_t> newConstr;
          std::set_difference(v.begin(), v.end(), active.begin(), active.end(),
                              std::inserter(newConstr, newConstr.begin()));

          if (newConstr.size()) // xeq in not feasible
            {
              Eigen::VectorXd xtemp=x0;
              std::set<size_t> va=Box_BestFeasible(x, xeq, box, 
                                                   newConstr, xtemp);              
              x=xtemp;
              std::copy(va.begin(), va.end(), 
                        std::inserter(active, active.begin()));
            }
          else
            {
              x=xeq;
              std::set<size_t> L=negativeIndices(lambda);
              if(L.size()==0)
                {
                  xout=x;
                  return 0;
                }
              else
                {
                  
                  std::set<size_t>::iterator f=active.begin();
                  for(size_t i=0; i<*L.begin(); ++i)
                    ++f;
                  active.erase(f);
                }
            }
      }
  }


}
