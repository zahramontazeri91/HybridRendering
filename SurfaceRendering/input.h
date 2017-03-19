#ifndef INPUT
#define INPUT

#include <Eigen/Dense>

using namespace Eigen;

VectorXd input_columns();
VectorXd input_rows();
MatrixXd input_pattern();
MatrixXd get_first_half_patch();
MatrixXd get_last_half_patch();
void get_aligned_masks();
int get_patch_number();


#endif
