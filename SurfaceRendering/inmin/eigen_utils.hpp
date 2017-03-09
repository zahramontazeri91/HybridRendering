/*
   (C) 2011 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

   \file eigen_utils.hpp

   Utilities for using the Eigen library
*/   
#ifndef _INMIN_EIGEN_UTILS_HPP__
#define _INMIN_EIGEN_UTILS_HPP__

#include <limits>

#include <Eigen/Core>
#include <Eigen/SVD>

namespace inmin {

  /** Compute the Pseudo-Inverse of a matrix 

      \param Ap The result is stored here

      \return The rank
   */
  inline size_t PseudoInv(const Eigen::MatrixXd &A,
			  Eigen::MatrixXd &Ap)
  {
    using namespace Eigen;

    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);


    JacobiSVD<MatrixXd>::SingularValuesType s=svd.singularValues();
    JacobiSVD<MatrixXd>::SingularValuesType inv=VectorXd::Zero(s.size());

    size_t rank;
    for (rank=0; rank< (size_t)A.cols(); ++rank) 
    {
      if (s[rank] > std::numeric_limits<double>::epsilon() *  s[0])
      {
	inv[rank]=1.0/s[rank];
      }
      else
      {
	break;
      }
    }
    Ap= svd.matrixV() *inv.asDiagonal()* svd.matrixU().transpose();
    return rank;
  }
}

#endif
