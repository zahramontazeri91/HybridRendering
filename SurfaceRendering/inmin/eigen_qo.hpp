/*
   (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>

   LICENSE: GNU General Public License V2  

   \file eigen_qo.hpp

   Quadratic optimisation using the Eigen library
*/   

#include <set>
#include <Eigen/Core>

using namespace Eigen;

/** Definition of constraints
   
   \note This should be moved to a public header 
 */
typedef struct 
{
  
} inmin_constraint;



namespace inmin {

  /*  Compute the basic constrained Qudratic Optimisation problem

      Constraints are defined by A^t x -b =0.
      
      \param H is the quadratic definition matrix
      \param g is the "gradient" vector
      \param A is the matrix part of the constraints
      \param b is the vector part of the constraints
   */
  int BQO_raw(const Eigen::MatrixXd &H,
               const Eigen::VectorXd &g,
               const Eigen::MatrixXd &A,
               const Eigen::VectorXd &b,
               Eigen::VectorXd &xout,
              Eigen::VectorXd &lambda);

  /** Compute the set of box constraints that is broken by a point
   */
  std::set<size_t> Box_broken(const Eigen::VectorXd &x,
                              const double *box);

  /** Compute the closest feasible point to xeq and constraints on
      which this point lies.
   */
  std::set<size_t> Box_BestFeasible(const Eigen::VectorXd &x0,
                                    const Eigen::VectorXd &xeq,
                                    const double *box,
                                    const std::set<size_t> &newConstr,
                                    Eigen::VectorXd &xout);

  /** Make matrix A and vector b defining equality constraints

      \param n Number of parameters
   */
  void Box_MakeConstrA(const std::set<size_t> &active,
                       const double *box,
                       size_t n,
                       Eigen::MatrixXd &A,
                       Eigen::VectorXd &b);
  

  /** Return indices of negative elements of x
   */
  std::set<size_t> negativeIndices(const Eigen::VectorXd &x);
  

  /* General quadratic optimisation with box constraints
     
     \param H is the quadratic definition matrix

     \param g is the "gradient" vector     

     \param box Defines the constraint box (see inmin_lm_in) for
     definition

     \param x0 is the initial guess 

   */
  int GQO_box(const Eigen::MatrixXd &H,
              const Eigen::VectorXd &g,
              const Eigen::VectorXd &x0,
              const double *box,
              Eigen::VectorXd &xout);

}
