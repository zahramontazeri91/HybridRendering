#ifndef INPUT
#define INPUT

#include <Eigen/Dense>

using namespace Eigen;

VectorXd input_columns();
VectorXd input_rows();
MatrixXd input_pattern();
void get_aligned_masks();

#endif
